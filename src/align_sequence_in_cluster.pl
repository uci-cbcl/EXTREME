#!/usr/bin/perl -w
use FindBin qw($Bin);
use lib "$Bin";
use gene_utils;

die "usage: perl $0 foo.sim foo.motif.list foo.cluster" if ($#ARGV<2);

@motifList = get_file_data($ARGV[1]); chomp @motifList;
@id2index = ();
for ($i=0; $i<=$#motifList; $i++) {
    @a = split(" ", $motifList[$i]);
    $id2index{$a[0]} = $i;
}

# similarity matrix
@all = get_file_data($ARGV[0]); chomp @all;
%simHash=();
foreach (@all) {
    @a = split(" ");
    $simHash{"$a[0]\t$a[1]"} = {'pos'=>$a[4],'dir'=>$a[3], 'score'=>$a[2]};
}

# cluster
@all = get_file_data($ARGV[2]); chomp @all;
for ($i=0; $i<=$#all; $i++) {
    next if (not $all[$i]=~/^>/);
    print $all[$i],"\n";
    
    @tmpList = ();
    while ($all[++$i]=~/[a-z]/i) {
	@a = split(" ", $all[$i]);
	push @tmpList, {'seq'=>$a[0], 'id'=>$a[0], 'ref'=>$all[$i], 'space'=>0, 'dir'=>'+', 'seqLen'=>length($a[0])};
    }


    
    for ($j=1; $j<=$#tmpList; $j++) {
	$idB = ${$tmpList[$j]}{'id'};
	die "Error $idB not exists" if (not exists $id2index{$idB});
	$indexB = $id2index{$idB};

	# find matched idA to align to
	for ($k=$j-1; $k>=0; $k--) {
	    $idA = ${tmpList[$k]}{'id'};
	    die "Error $idA not exists" if (not exists $id2index{$idA});
	    $indexA = $id2index{$idA};
	    if ($indexA < $indexB) {
		next if (not exists $simHash{"$indexA\t$indexB"});
	    } else {
		next if (not exists $simHash{"$indexB\t$indexA"});
	    }
	    last;  # break of found one
	}

	
	if ($k<0) {
	    foreach (@tmpList) {
		print $_->{'id'},"\t", $id2index{$_->{'id'}},"\n";
	    }
	}

	die "Error: $idA not found any sim" if ($k<0);
	
	
	#print "### align two: \n";
	#print ${$tmpList[$j]}{'ref'},"\n";
	#print ${$tmpList[$k]}{'ref'},"\n";



	$alignDirB = '';
	$alignPosB = '';
	
	if ($indexA < $indexB) {
	    ($pos, $dir) = (${$simHash{"$indexA\t$indexB"}}{'pos'}, ${$simHash{"$indexA\t$indexB"}}{'dir'});
	    
	    #print "sim: $pos $dir\n";

	    $alignDirB = $dir;
	    $alignPosB = $pos;

	    if (${$tmpList[$k]}{'dir'} eq '-') {
		$alignDirB = ($dir eq '+' ? '-' : '+');
		$alignPosB =  ${$tmpList[$k]}{'seqLen'} - $pos - ${$tmpList[$j]}{'seqLen'};
	    }

	} else {
	    ($pos, $dir) = (${$simHash{"$indexB\t$indexA"}}{'pos'}, ${$simHash{"$indexB\t$indexA"}}{'dir'});

	    #print "sim: $pos $dir\n";

	    $alignDirB = $dir;
	    $alignPosB = $pos;
	    
	    if (${$tmpList[$k]}{'dir'} eq '+') {
		if ($dir eq '+') {
		    $alignPosB = - $pos;
		} else {
		    $alignPosB = - (${$tmpList[$j]}{'seqLen'} - $pos - ${$tmpList[$k]}{'seqLen'});
		}
	    } else {
		if ($dir eq '+') {
		    $alignDirB = '-';
		    $alignPosB = - (${$tmpList[$j]}{'seqLen'} - $pos - ${$tmpList[$k]}{'seqLen'});
		} else {  # already inversed
		    $alignDirB = '+';
		    $alignPosB = - $pos;
		}
	    }
	}
	


	${$tmpList[$j]}{'seq'} = reverse_comp(${$tmpList[$j]}{'seq'}) if ($alignDirB eq '-');
	${$tmpList[$j]}{'dir'} = $alignDirB;
	
	$alignPosB += ${$tmpList[$k]}{'space'};

	#print "### final align $alignDirB  $alignPosB\n";

	if ($alignPosB > 0) {
	    ${$tmpList[$j]}{'space'} += $alignPosB;
	} elsif ($alignPosB < 0) { # add space before all previous motifs
	    for ($k=0; $k<$j; $k++) {
		${$tmpList[$k]}{'space'} -= $alignPosB;
	    }
	}
	
    }
    
    foreach (@tmpList) {
	print '.' x ${$_}{'space'},  ${$_}{'seq'}, "\t", ${$_}{'ref'},"\t",${$_}{'dir'},"\n";
    }
    
    print "\n";
}

    
