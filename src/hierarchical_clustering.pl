#!/usr/bin/perl -w
use FindBin qw($Bin);
use lib "$Bin";
use gene_utils;


die "usage: perl NAME dist_matrix threshold_dist default_missing_distance" if ($#ARGV<2);

$dist_matrix_file = $ARGV[0];
$threshold = $ARGV[1];
$missingPenalty = $ARGV[2];

## read motif pair simalirity file
# first pass to read id first:
%tmpHash=();
open IN, $dist_matrix_file or die "cannot open $dist_matrix_file";
while (<IN>) {
    @a = split(" ");
    $tmpHash{$a[0]}=1;
    $tmpHash{$a[1]}=1;
}
close IN;
@idList = keys %tmpHash; 

# sort the id index, not necessary
@idList = sort {$a<=>$b} @idList; 
%id2index=();
for ($i=0; $i<@idList; $i++) {
    $id2index{$idList[$i]} = $i;
}

## read the sim file in the second pass:, read dist matrix
%distMatrix=(); @pairList=();
open IN, $dist_matrix_file or die "cannot open $dist_matrix_file";
while (<IN>) {
    @a = split(" ");
    die "Error $a[0] , $a[1] not exists" if ((not exists $id2index{$a[0]}) or (not exists $id2index{$a[1]}));
    $idA = $id2index{$a[0]};
    $idB = $id2index{$a[1]};
    ($idA,$idB) = ($idB,$idA) if ($idA > $idB);  # always make sure $idA < $idB
    $distMatrix{"$idA\t$idB"} = $a[2];
    push @pairList, {'A'=>$idA, 'B'=>$idB, 'score'=>$a[2]} if ($a[2] <= $threshold);
}
close IN;

# so from now, everything is referenced to index, rather than the original ID

# to store merged ids
@clusterList = (); $#clusterList = $#idList;
for ($i=0; $i<@clusterList; $i++) {
    push @{$clusterList[$i]}, $i;
}


# iterative sort
while (scalar(@pairList) > 0) {

    $t = merge_nearest_pair();
    last if ($t < 1);
}




# output list
@outList = ();
foreach $c (@clusterList) {
    next if (scalar(@{$c}) < 1);
    @a = ();
    foreach (@{$c}) {
	push @a, $idList[$_];
    }
    push @outList, {'list'=>join("\n",@a),'num'=>scalar(@a)};
}
@outList = sort {$b->{'num'}<=>$a->{'num'}} @outList;
$k = 0;
foreach (@outList) {
    $k++;
    print ">cluster$k\t", $_->{'num'},"\n";
    print $_->{'list'},"\n\n";
}




sub merge_nearest_pair {

    # find nearest pair, smallest score
    my $bestIndex=0; my $bestScore=$threshold+1;
    for (my $i=0; $i<@pairList; $i++) {
	if ($pairList[$i]->{'score'} < $bestScore) {
	    $bestScore = $pairList[$i]->{'score'};
	    $bestIndex = $i;
	}
    }
    
    return 0 if ($bestScore > $threshold);
    
    my $idA = $pairList[$bestIndex]->{'A'};
    my $idB = $pairList[$bestIndex]->{'B'};
    splice @pairList, $bestIndex, 1;
    
    # store previous number of elements in the two clusters
    $numA = scalar(@{$clusterList[$idA]});
    $numB = scalar(@{$clusterList[$idB]});
    
    # merger all cluster B -> cluster A
    push @{$clusterList[$idA]}, @{$clusterList[$idB]};
    @{$clusterList[$idB]} = (); 
    
    # update pairwise distance matrix
    # update 1. update idB -> idA

    # step 1: store previous other nodes shared with idA, store score
    my %tmpHashA=();
    for (my $j=0; $j<@pairList; $j++) {
	$tA = $pairList[$j]->{'A'};
	$tB = $pairList[$j]->{'B'};
	if ($tA == $idA) {
	    $tmpHashA{$tB} = {'pairIndex'=>$j, 'updated'=>0};
	} elsif ($tB == $idA) {
	    $tmpHashA{$tA} = {'pairIndex'=>$j, 'updated'=>0};
	}
    }
    
    
    my @toRemovePairList = (); my $score=''; my $i='';
    for (my $j=0; $j<@pairList; $j++) {
	$tA = $pairList[$j]->{'A'};
	$tB = $pairList[$j]->{'B'};
	next if ($tA != $idB and $tB != $idB);
	($tA, $tB) = ($tB, $tA) if ($tA == $idB);  # always in the format of tB==tdB
	
	if (exists $tmpHashA{$tA}) {
	    # update the pair store A, and remove current j
	    $i = $tmpHashA{$tA}->{'pairIndex'};
	    $score = ($pairList[$i]->{'score'} * $numA + $pairList[$j]->{'score'} * $numB)/($numA+$numB);
	    if ($score > $threshold) {
		push @toRemovePairList, $i;
	    } else {
		$pairList[$i]->{'score'} = $score;
	    }
	    $tmpHashA{$tA}->{'updated'} = 1;
	    push @toRemovePairList, $j;
	} else {
	    # update, just change idB to idA
	    $score = ($pairList[$j]->{'score'}*$numB + $missingPenalty*$numA)/($numA+$numB);
	    if ($score > $threshold) {
		push @toRemovePairList, $j;
	    } else {
		# keep the node, turn to right order
		if ($tA <= $idA) {
		    $pairList[$j]->{'B'} = $idA;
		} else {
		    $pairList[$j]->{'A'} = $idA;
		    $pairList[$j]->{'B'} = $tA;
		}
		$pairList[$j]->{'score'} = $score;
	    }
	}
    }
    
    # step 3. process those pairs with only idA, not no coexisting with idB
    my @t = keys %tmpHashA;
    foreach (@t) {
	next if ($tmpHashA{$_}->{'updated'} == 1);

	$i = $tmpHashA{$_}->{'pairIndex'};
	$score = ($pairList[$i]->{'score'}*$numA + $missingPenalty*$numB)/($numA+$numB);
	if ($score > $threshold) {
	    push @toRemovePairList, $i;
	} else {
	    $pairList[$i]->{'score'} = $score;
	}
    }

    @toRemovePairList = sort {$b<=>$a} @toRemovePairList;
    foreach (@toRemovePairList) {
	splice @pairList, $_, 1;
    }

    return 1;
}


# unweighted average distance
sub calculate_distance_between_two_groups_average {
    my ($idListARef, $idListBRef) = @_;
    my $dist = 0;
    my $n = 0;
    my ($a, $b) = (); my $s='';
    foreach (@$idListARef) {
	$a = $_;
	foreach  (@$idListBRef) {
	    $b = $_;
	    next if ($b == $a);
	    $n++;

	    if ($a <= $b) {
		$s = (exists $distMatrix{"$a\t$b"} ? $distMatrix{"$a\t$b"} : $missingPenalty);
	    } else {
		$s = (exists $distMatrix{"$b\t$a"} ? $distMatrix{"$b\t$a"} : $missingPenalty);
	    }

	    $dist += $s;
	}
    }
    return $dist/$n;
}

