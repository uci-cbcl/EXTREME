
%genetic_code = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
		);



# get_file_data
# A subroutine to get data from a file given its filename

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}



# write a string list to a file
sub write_list_to_file {
    my ($filename, @a) = @_;

    use strict;
    use warnings;
    
    unless( open(tmp_out, ">$filename") ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    chomp @a;
    foreach (@a) {
	print tmp_out $_,"\n";
    }
    close tmp_out;
}




# Translate IUB ambiguity codes to regular expressions 
# IUB_to_regexp
#
# A subroutine that, given a sequence with IUB ambiguity codes,
# outputs a translation with IUB codes changed to regular expressions
#
# These are the IUB ambiguity codes
# (Eur. J. Biochem. 150: 1-5, 1985):
# R = G or A
# Y = C or T
# M = A or C
# K = G or T
# S = G or C
# W = A or T
# B = not A (C or G or T)
# D = not C (A or G or T)
# H = not G (A or C or T)
# V = not T (A or C or G)
# N = A or C or G or T 

sub IUB_to_regexp {

    my($iub) = @_;

    $iub = uc $iub;
    
    use strict;
    use warnings;

    my $regular_expression = '';

    my %iub2character_class = (
    
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
	N => '[ACGT]',
	n => '[ACGT]',		       
    );


    # Translate each character in the iub sequence
    for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
        $regular_expression
          .= $iub2character_class{substr($iub, $i, 1)};
    }

    return $regular_expression;
}


sub IUB_to_regexp_ext {

    my($iub) = @_;
    
    $iub = uc $iub;

    use strict;
    use warnings;

    my $regular_expression = '';

    my %iub2character_class = (
    
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
        N => '[ACGT]',
        n => '[ACGT]',
    );

    # Translate each character in the iub sequence
    for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
        $regular_expression
          .= "$iub2character_class{substr($iub, $i, 1)}-*";
    }
    $regular_expression=~s/-\*$//;

    return $regular_expression;
}




#
# Find locations of a match of a regular expression in a string
#
# 
# return an array of positions where the regular expression
#  appears in the string
#
# position start from 0
sub match_positions {

    my($regexp, $sequence) = @_;

    use strict;
    use warnings;

    #
    # Declare variables
    #

    my @positions = ( );

    #
    # Determine positions of regular expression matches
    #
    
    while ( $sequence =~ /$regexp/ig ) {

        push ( @positions, pos($sequence) - length($&) );
    }

    return @positions;
}


sub match_positions_case_sensitive {

    my($regexp, $sequence) = @_;

    use strict;
    use warnings;

    #
    # Declare variables
    #

    my @positions = ( );

    #
    # Determine positions of regular expression matches
    #
    
    while ( $sequence =~ /$regexp/g ) {

        push ( @positions, pos($sequence) - length($&) );
    }

    return @positions;
}


# return sorted index in ascending order
sub index_sort_ascend {
    return reverse(index_sort(@_));
}

# return sorted index in descending order
sub index_sort_descend {
    return index_sort(@_);
}

# sort function 
# return sorted index in descending order
sub index_sort {

    my @list = @_;

    use strict;
    use warnings;

    my @index=( );
    my @sorted_index=( );
    my $i='';

    for ($i=0;$i<=$#list;$i++) {
	$index[$i]=$i;
    }
    
    @sorted_index=sort {$list[$b] <=> $list[$a]} @index;
    
    return @sorted_index;
}


# find min index;
sub index_min { 
    use strict;
    use warnings;
    my @list = @_;
    my @sorted_index = index_sort_descend(@list);
    return pop @sorted_index;
}


# find max index;
sub index_max {
    use strict;
    use warnings;
    my @list = @_;
    my @sorted_index = index_sort_descend(@list);
    return shift @sorted_index;
}

# find max value;
sub max {
    my @list = @_;
    use strict;
    use warnings;

    my $m =  $list[0];
    for (my $i=1;$i<=$#list; $i++) {
	$m = $list[$i] if ($list[$i]>$m);
    }
    return $m;
}


# find min value;
sub min {
    my @list = @_;
    use strict;
    use warnings;

    my $m =  $list[0];
    for (my $i=1;$i<=$#list; $i++) {
	$m = $list[$i] if ($list[$i]<$m);
    }
    return $m;
}








# compare two list, return index of equal elements
# pass by reference
# return: shared_element   index_from_array_a index_from_array_b
sub array_comp {
    my ($a, $b) = @_;
    
    use strict;
    use warnings;
    
    my $i='';
    my @list_a = @$a;
    my @list_b = @$b;
    my @return_list = ( );
    my @idx = ( );

    
    for($i=0; $i<=$#list_a; $i++) {
	@idx = index_array(\@list_b, $list_a[$i]);
	if (@idx) {
	    foreach (@idx) {
		push @return_list, "$list_a[$i] \t $i \t $_";
	    }
	}
    }
    return @return_list;
}



# return index of an element in an array
# pass array by reference
sub index_array {
    my ($p, $element) = @_;
    
    use strict;
    use warnings;

    my @list = ( );
    my $i = '';

    for( my $i=0; $i<=$#$p; $i++) {
	if ($$p[$i] eq $element) {
	    push @list, $i;
	}
    }
    
    return @list;
}


sub reverse_comp {
    my ($seq) = @_;

    $seq = reverse $seq;
    $seq=~tr/ATCGatcg/TAGCtagc/;
    $seq=~tr/RYKM/YRMK/;
    $seq=~tr/BDHV/VHDB/;
    
    return $seq;
}




# get aligned ace sequences from a ACE file
# input: ace filename, motif number
sub get_aligned_ace_seqs {
    use strict;
    use warnings;

    my ($ace_filename, $motif_index) = @_;
    
    my @all_data = get_file_data($ace_filename);
    chomp @all_data;

    my @seqs_vec = ();

    for (my $i=0; $i<=$#all_data; $i++) {
	if (not $all_data[$i] =~ /^Motif/) {next; } 
	my @a = split(" ", $all_data[$i]);
	my $b = pop @a;
	if ($b != $motif_index) {next;}
    
	while (not $all_data[++$i] =~ /\*/) {
	    @a=split(" ",$all_data[$i]);
	    push @seqs_vec, $a[0];
	}
	
	$i=$#all_data+1;
    }
    
    return @seqs_vec;
}




# get aligned ace sequences from a ACE file
# input: ace filename, motif number
# output: fasta format
sub get_aligned_ace_seqs_4_fasta {
    use strict;
    use warnings;

    my ($ace_filename, $motif_index) = @_;
    
    my @all_data = get_file_data($ace_filename);
    chomp @all_data;

    my @seqs_vec = ();

    for (my $i=0; $i<=$#all_data; $i++) {
	if (not $all_data[$i] =~ /^Motif/) {next; } 
	my @a = split(" ", $all_data[$i]);
	my $b = pop @a;
	if ($b != $motif_index) {next;}
    
	while (not $all_data[++$i] =~ /\*/) {
	    @a=split(" ",$all_data[$i]);
	    my $seq = shift @a;
	    my $otherinfo = join("\t",@a);
	    push @seqs_vec, ("> $ace_filename\t$motif_index\t$otherinfo", $seq);
	}
	
	$i=$#all_data+1;
    }
    
    return @seqs_vec;
}

sub calculate_entropy_given_nucleotide_string_normalized {
    my ($seq) = @_;
    return calculate_entropy_given_nucleotide_string($seq)/log(4)*log(2);
}

# calculate entropy of a nucleotide sequence ACGT only, given seq string
sub calculate_entropy_given_nucleotide_string {
    my ($seq) = @_;
    use strict;
    use warnings;
    my %tmpHash=(); 
    my @a = split("", uc($seq)); 
    foreach (@a) {
	$tmpHash{$_}++;
    }
    my @baseList = ('A','C','G','T');
    my @count = (); $#count=3;

    for (my $i=0; $i<4; $i++) {
	if (exists $tmpHash{$baseList[$i]}) {
	    $count[$i] = $tmpHash{$baseList[$i]};
	} else {
	    $count[$i] = 0;
	}
    }
    # add pseudo count
    for (my $i=0; $i<=$#count; $i++) {
	$count[$i] += 0.01;
    }
    return calc_entropy(@count);
}


# get entropy of 2mer distribution in a given sequence
sub calculate_entropy_given_nucleotide_string_2mer_normalized {
    my ($seq) = @_;
    return calculate_entropy_given_nucleotide_string_2mer($seq)/log(16)*log(2);
}

# get entropy of 2mer distribution in a given sequence
sub calculate_entropy_given_nucleotide_string_2mer {
    my ($seq) = @_;
    use strict;
    use warnings;
    $seq = uc $seq;
    my %tmpHash=(); my $L = length($seq)-2;
    for (my $i=0; $i<=$L; $i++) {
	$tmpHash{substr($seq,$i,2)}++;
    }
    
    my @baseList = ('A','C','G','T');
    my @count = (); $#count=15;
    my ($b1, $b2) = ();
    my $n = 0;
    foreach $b1 (@baseList) {
	foreach $b2 (@baseList) {
	    if (exists $tmpHash{$b1.$b2}) {
		$count[$n] = $tmpHash{$b1.$b2};
	    } else {
		$count[$n] = 0;
	    }
	    $n++;
	}
    }
    # add pseudo count
    for (my $i=0; $i<=$#count; $i++) {
	$count[$i] += 0.01;
    }
    
    return calc_entropy(@count);
}


# calculate entropy of a prob distribution, log2 based
sub calc_entropy_for_prob_row_log2 {
    use strict;
    use warnings;
    my @p = @_;
    # normalize
    my $sum = 0;
    foreach (@p) {
	$_ += 0.000001;  # add a small eps term
	$sum += $_;
    }
    foreach (@p) {
	$_ = $_/$sum;
    }
    
    my $log2 = log(2);
    my $H = 0;
    foreach (@p) {
	$H -= $_*log($_)/$log2;
    }
    return $H;
}


# calculate entropy of a prob distribution
sub calc_entropy {
    use strict;
    use warnings;
    my @p = @_;
    # normalize
    my $sum = 0;
    foreach (@p) {
	$_ += 0.00001;  # add a small eps term
	$sum += $_;
    }
    foreach (@p) {
	$_ = $_/$sum;
    }
    
    my $log2 = log(2);
    my $H = 0;
    foreach (@p) {
	$H -= $_*log($_)/$log2;
    }
    return $H;
}





# compute information of a prob dist
sub calc_info {
    use strict;
    use warnings;

    my @p = @_;
    my $log2 = log(2);
    my $info = 0;
    
    # normalize
    my $sum = 0;
    for (my $i=0; $i<=$#p; $i++) {
	$sum += $p[$i];
    }
    for (my $i=0; $i<=$#p; $i++) {
	$p[$i]= ($p[$i]+0.0001)/($sum+0.0004);
    }
    

    for (my $i=0; $i<=$#p; $i++) {
	$info += - $p[$i]*log($p[$i])/$log2;
    }
    return $info;
}



# get 2-information for DNA sequences
sub get_order_4_seqs {
    use strict;
    use warnings;
    
    my @order = get_info_4_seqs(@_);
    for (my $i=0; $i<=$#order; $i++) {
	$order[$i] = 2 - $order[$i];
    }
    return @order;
}


# compute information for DNA sequences
sub get_info_4_seqs {
    use strict;
    use warnings;

    my @seqs = @_;

    my @stat = ();
    my @info = ();
    my $L = '';

    for (my $j=0;$j<length($seqs[0]);$j++) {
	my $c='';
	for (my $k=0;$k<=$#seqs;$k++) {
	    $c.=substr($seqs[$k],$j,1);
	}
	$L = ($c=~tr/actgACTG//);
	$stat[$j][0]=($c=~tr/aA// + 0.1)/($L + 0.4);
	$stat[$j][1]=($c=~tr/tT// + 0.1)/($L + 0.4);
	$stat[$j][2]=($c=~tr/cC// + 0.1)/($L + 0.4);
	$stat[$j][3]=($c=~tr/gG// + 0.1)/($L + 0.4);

	$info[$j] = calc_info(@{$stat[$j]});
    }

    return @info;
}

# get raw count of bases in aligned sequences
# row i represents counts of A,C,G,T at position i respectively
sub get_count_matrix_given_seq_list {
    use strict;
    use warnings;

    my @seqs = @_;
    my @stat = ();

    for (my $j=0;$j<length($seqs[0]);$j++) {
	my $c='';
	for (my $k=0;$k<=$#seqs;$k++) {
	    $c.=substr($seqs[$k],$j,1);
	}
	$stat[$j][0]=($c=~tr/[aA]//);
	$stat[$j][1]=($c=~tr/[cC]//);
	$stat[$j][2]=($c=~tr/[gG]//);
	$stat[$j][3]=($c=~tr/[tT]//);
    }
    return @stat;
}


# get raw count of bases in aligned sequences
# row i represents counts of A,C,G,T at position i respectively
# but return probability distribution
# input sequences reference and pseudo counts
sub get_weight_matrix_given_seq_list_return_prob_matrix {
    use strict;
    use warnings;

    my ($seqRef, $pseudoCount) = @_;
    my @count = get_count_matrix_given_seq_list(@$seqRef);
    
    my $s=''; my $r='';
    foreach $r (@count) {
	$s = 0;
	foreach (@{$r}) {
	    $_ += $pseudoCount;
	    $s += $_;
	}
	foreach (@{$r}) {
	    $_ = $_/$s;
	}
    }
    return @count;
}


# get raw count of bases in aligned sequences
# row i represents counts of A,C,G,T at position i respectively
# but return probability distribution
# input sequences reference and pseudo counts
# output string lines, with 4 prob, consensus and info
sub get_weight_matrix_given_seq_list_return_prob_and_consensus_info_string {
    use strict;
    use warnings;

    my ($seqRef, $pseudoCount) = @_;
    my @count = get_count_matrix_given_seq_list(@$seqRef);
    
    my $s=''; my $r='';
    my @outList = (); my $c=''; my $info='';
    foreach $r (@count) {
	$s = 0;
	foreach (@{$r}) {
	    $_ += $pseudoCount;
	    $s += $_;
	}
	foreach (@{$r}) {
	    $_ = $_/$s;
	}

	$c = get_consensus_for_one_row_prob(@{$r});
	$info = 2 - calc_entropy_for_prob_row_log2(@{$r});
	$info = sprintf('%.2f', $info);

	foreach (@$r) {
	    $_ = sprintf('%.4f',$_);
	}
	push @outList, join("\t", (@$r, $c, $info));
    }
    return @outList;
}



sub get_weight_matrix_given_count_matrix_return_prob_and_consensus_info_string {
    use strict;
    use warnings;

    my @count = @_;
    my $pseudoCount = 0.01;
    
    my $s=''; my $r='';
    my @outList = (); my $c=''; my $info='';
    foreach $r (@count) {
	$s = 0;
	foreach (@{$r}) {
	    $_ += $pseudoCount;
	    $s += $_;
	}
	foreach (@{$r}) {
	    $_ = $_/$s;
	}

	$c = get_consensus_for_one_row_prob(@{$r});
	$info = 2 - calc_entropy_for_prob_row_log2(@{$r});
	$info = sprintf('%.2f', $info);

	foreach (@$r) {
	    $_ = sprintf('%.4f',$_);
	}
	push @outList, join("\t", (@$r, $c, $info));
    }
    return @outList;
}








# get weight matrix from aligned sequences
sub get_weight_matrix {
    use strict;
    use warnings;

    my @seqs = @_;
    my @stat = ();

    my @bases=('A','C','G','T');

    for (my $j=0;$j<length($seqs[0]);$j++) {
	my $c='';
	for (my $k=0;$k<=$#seqs;$k++) {
	    $c.=substr($seqs[$k],$j,1);
	}
	$stat[0][$j]=($c=~tr/[aA]//);
	$stat[1][$j]=($c=~tr/[cC]//);
	$stat[2][$j]=($c=~tr/[gG]//);
	$stat[3][$j]=($c=~tr/[tT]//);
    }
    
    my @output = ();
    for (my $i=0; $i<=$#stat; $i++) {
	my $b = $bases[$i];
	push @output,  $b."\t". join("\t",@{$stat[$i]}) ;
    }

    return @output;
}



# get weight matrix from aligned sequences
sub get_weight_matrix_and_info {
    use strict;
    use warnings;

    my @seqs = @_;
    my @stat = ();

    my @bases=('A','C','G','T');

    for (my $j=0;$j<length($seqs[0]);$j++) {
	my $c='';
	for (my $k=0;$k<=$#seqs;$k++) {
	    $c.=substr($seqs[$k],$j,1);
	}
	$stat[$j][0]=($c=~tr/[aA]//);
	$stat[$j][1]=($c=~tr/[cC]//);
	$stat[$j][2]=($c=~tr/[gG]//);
	$stat[$j][3]=($c=~tr/[tT]//);
    }

    my @output = ();
    for (my $i=0; $i<=$#stat; $i++) {
	my $info = 2 - calc_info(@{$stat[$i]});
	my $index = index_max(@{$stat[$i]});
	my $b = $bases[$index];
	push @output,  $b."\t". join("\t",@{$stat[$i]}) ."\t".$info;
    }

    return @output;
}


# given prob of A,C,G,T respectively
sub get_consensus_for_one_row_prob {
    use strict;
    use warnings;
    my ($pA, $pC, $pG, $pT) = @_;
    my $s = $pA+$pC+$pG+$pT;
    $pA = $pA/$s;
    $pC = $pC/$s;
    $pG = $pG/$s;
    $pT = $pT/$s;

    if ($pA>=0.7) { 
	return 'A'; 
    } elsif ($pT>=0.7) { 
	return 'T'; 
    } elsif ($pC>=0.7) { 
	return 'C'; 
    } elsif ($pG>=0.7) { 
	return 'G'; 
    } elsif ($pA+$pT>0.8) {
	return 'W';
    } elsif ($pA+$pC>0.8) {
	return 'M';
    } elsif ($pA+$pG>0.8) {
	return 'R';
    } elsif ($pT+$pC>0.8) {
	return 'Y';
    } elsif ($pT+$pG>0.8) {
	return 'K';
    } elsif ($pC+$pG>0.8) {
	return 'S';
    } elsif ($pA+$pT>0.8) {
	return 'W';
    } elsif ($pA<0.1) {
	return 'B';
    } elsif ($pT<0.1) {
	return 'V';
    } elsif ($pC<0.1) {
	return 'D';
    } elsif ($pG<0.1) {
	return 'H';
    } else {
	return 'N';
    }
}



# get consensus sequence from weight matrix input
# weight matrix is in the format of (2-dim array)
#    A C G T
# p1 . . . .
# p2 . . . .
# p3 . . . .
# p4 . . . .
sub get_consensus_from_wm_row {
    use strict;
    use warnings;
    my @wm = @_; my ($pA,$pC,$pG,$pT)=(); my $consensus='';
    foreach (@wm) {
	($pA, $pC, $pG, $pT) = ($_->[0],$_->[1],$_->[2],$_->[3]);
	$consensus .= get_consensus_for_one_row_prob($pA,$pC,$pG,$pT);
    }
    return $consensus;
}

# get consensus sequence from weight matrix input
# weight matrix is in the format of 
# A num1 num2 ... 
# C ...  ...  ...
# G ...  ...  ...
# T ...  ...  ...
sub get_consensus_from_wm {
    
    use strict;
    use warnings;

    my @weight = @_;
    my @A=(); my @T=(); my @G=(); my @C=();

    foreach (@weight) {
	my @tmp = split(" ", $_);
	my $first = shift @tmp;
	
	if ($first eq 'A' or $first eq 'a') {
	    @A = @tmp; 
	} elsif ($first eq 'T' or $first eq 't') {
	    @T = @tmp;
	} elsif ($first eq 'C' or $first eq 'c') {
	    @C = @tmp;
	} elsif ($first eq 'G' or $first eq 'g') {
	    @G = @tmp;
	}
    }
    
    my $consensus = '';
    for (my $i=0; $i<=$#A; $i++) {
	my $total_num = $A[$i]+$T[$i]+$C[$i]+$G[$i];
	my $pA=$A[$i]/$total_num;
	my $pT=$T[$i]/$total_num;
	my $pC=$C[$i]/$total_num;
	my $pG=$G[$i]/$total_num;
	$consensus .= get_consensus_for_one_row_prob($pA,$pC,$pG,$pT);
    }
    
    #$consensus =~ s/^[BVDHN]+//;
    #$consensus =~ s/[BVDHN]+$//;
    
    return $consensus;
}




# get consensus sequence aligned sequences
# each aligned sequence is in a row
sub get_consensus_from_seqs {

    use strict;
    use warnings;

    my @seqs = @_;
    my @weight = get_weight_matrix(@seqs);
    my $consensus = get_consensus_from_wm(@weight);
    
    return $consensus;
}


# return the number of same chars between two strings
sub get_identity_between_strs {
    use strict;
    use warnings;
    
    my ($str1, $str2) = @_;
    my $n = 0;
    for (my $i=0; $i<length($str1); $i++) {
	if (substr($str1,$i,1) eq substr($str2,$i,1)) {
	    $n++;
	}
    }

    return $n;
}




# return log(n!), approximated using Sterling's Formula
sub log_factorial {
    my ($n) = @_;
    use strict;
    use warnings;

    my $z = 0;
    if ($n<1) {
	return 0;
    } elsif ($n<5) {
	while ($n>1) {$z += log($n--);}
    } else {
	$z = $n*log($n)-$n+log(2*3.14159*$n)/2;
    }
    return $z;
}


# return hypergeometric distribution
sub geometric_prob {
    my ($N,$n,$K,$k) = @_;
    use strict;
    use warnings;

    my $p = log_factorial($K)+log_factorial($N-$K)+log_factorial($N-$n)+log_factorial($n)-log_factorial($N)-log_factorial($K-$k)-log_factorial($k)-log_factorial($N-$K-$n+$k)-log_factorial($n-$k);
    
    return exp($p);
}

sub get_left_tail_geometric_logpvalue {
    my ($N,$n,$K,$k) = @_;
    use strict;
    use warnings;

    my $pvalue = get_left_tail_geometric_pvalue($N,$n,$K,$k);
    return -30 if ($pvalue <=0);
    return log10($pvalue);
}


# return pvalue of hypergeometric distribution
sub get_left_tail_geometric_pvalue {
    my ($N,$n,$K,$k) = @_;
    use strict;
    use warnings;
    
    my $pvalue=0;
    for (my $i=0;$i<=$k;$i++) {
	$pvalue += geometric_prob($N,$n,$K,$i);
    }
    
    return 0 if ($pvalue<0);
    return $pvalue;
    
    return -10  if ($pvalue<=0);
    return log10($pvalue);
}



sub get_geometric_logpvalue {
    my ($N,$n,$K,$k) = @_;
    use strict;
    use warnings;

    my $pvalue = get_geometric_pvalue($N,$n,$K,$k);
    return -30 if ($pvalue <=0);
    return log10($pvalue);
}


# return pvalue of hypergeometric distribution
sub get_geometric_pvalue {
    my ($N,$n,$K,$k) = @_;
    use strict;
    use warnings;
    
    my $max_k = $n;
    if ($K<$n) {$max_k=$K;}
  
    my $pvalue=0;
    for (my $i=$k;$i<=$max_k;$i++) {
	$pvalue += geometric_prob($N,$n,$K,$i);
    }
    
    return 0 if ($pvalue<0);
    return $pvalue;

    return -10  if ($pvalue<=0);
    return log10($pvalue);
}


sub log10 {
    use strict;
    use warnings;

    my $n = shift;
    return log($n)/log(10);
}


# generate repeated string 
sub repeat_chr {
    my ($c, $n) = @_;
    use strict;
    use warnings;

    my $output='';
    while ($n>0) { $output .= $c; $n--; }
    return $output;
}


# generate Kmer list
# given num and/or base list
sub generate_Kmer_list {
    use strict;
    use warnings;

    my $num_mer = '';
    my @bases = '';
    
    if (scalar @_ < 2) {
	($num_mer) = @_;
	@bases = ('A','C','G','T');
    } else {
	($num_mer, @bases) = @_;
    }

  
    my @motif_list = ();

    
    if ($num_mer>0) {
	@motif_list = @bases;
    }
    
    my @a = (); 
    my ($m,$b) = ();
    for (my $i=0; $i<$num_mer-1; $i++) {
    	@a = ();
	foreach $m (@motif_list) {
	    foreach $b (@bases) {
		push @a, $m . $b;
	    }
	}
	@motif_list = @a;
    }

    return @motif_list;
}


sub sum_given_ref {
    my ($r) = @_;
    use strict;
    use warnings;
    
    my $s = 0;
    foreach (@{$r}) {
	$s += $_;
    }
    return $s;
}


sub std_given_ref {
    my ($a) = @_;
    use strict;
    use warnings;
    
    my $mu = mean_given_ref($a);
    
    my $s = 0;
    foreach (@{$a}) {
	$s += ($_-$mu)*($_-$mu);
    }

    my $n = scalar @{$a};
    return sqrt($s/$n);
}



sub mean_given_ref {
    my ($a) = @_;
    use strict;
    use warnings;
    my $n = scalar @{$a};
    if ($n<=0) {return 0;}
    else {
	return sum_given_ref($a)/$n;
    }
}

sub sum {
    my (@a) = @_;
    use strict;
    use warnings;
    
    my $s = 0;
    foreach (@a) {
	$s += $_;
    }
    return $s;
}


sub std {
    my (@a) = @_;
    use strict;
    use warnings;
    
    my $mu = mean(@a);
    
    my $s = 0;
    foreach (@a) {
	$s += ($_-$mu)*($_-$mu);
    }

    my $n = scalar @a;
    return sqrt($s/$n);
}



sub mean {
    my (@a) = @_;
    use strict;
    use warnings;
    my $n = scalar @a;
    if ($n<=0) {return 0;}
    else {
	return sum(@a)/$n;
    }
}

sub median {
    my (@t) = @_;
    use strict;
    use warnings;
    @t = sort {$a<=>$b} @t;
    my $n = $#t;
    if (int($n/2) == $n/2) {
	return $t[$n/2];
    } else {
	return ($t[int($n/2)]+$t[int($n/2)+1])/2;
    }
}


sub median_given_ref {
    my ($reft) = @_;
    use strict;
    use warnings;
    my @t=();
    @t = sort {$a<=>$b} @{$reft};
    my $n = $#t;
    if (int($n/2) == $n/2) {
	return $t[$n/2];
    } else {
	return ($t[int($n/2)]+$t[int($n/2)+1])/2;
    }
}








# return the module of $x/$y;
sub mod {
    my ($x, $y) = @_;
    use strict;
    use warnings;

    my $i = sprintf("%d", $x/$y);
    return $x-$i*$y;
}


sub toLowerCase {
    my ($x) = @_;
    use strict;
    use warnings;

    if ($x) {
	$x =~ tr/[A-Z]/[a-z]/;
    }
    return $x;
}


sub toUpperCase {
    my ($x) = @_;
    use strict;
    use warnings;
    
    $x =~ tr/[a-z]/[A-Z]/ if ($x);
    return $x;
}



## get mann_whitney rank_sum statistics zscore
# $g: set_num
# $N: total_num_of_ranked_list
# $rs: rank_sum_score
sub mann_whitney_zscore {
    my ($g, $N, $rs) = @_;

    my $mu = $g*($N+1)/2;
    my $sig = sqrt($g*($N-$g)*(1+$N)/12);
    my $z = ($sig>0 ? ($mu-$rs)/$sig : 0);
 
    return $z;
}



# calculating Komogorov-Smirnov p value
# $D is the max ES
# $Ne is the effective number $Ne=N1*N2/(N1+N2)
sub kstwo {
    my ($D, $Ne) = @_;
    
    return probks( (sqrt($Ne) + 0.12 + 0.11/sqrt($Ne))*$D );
}


## calculating QKS in Komogovor-Smirnov statistics, ref NR in C
sub probks {
    use strict;
    use warnings;

    my ($alam) = @_;
    
    my $EPS1 = 0.001;
    my $EPS2 = 1.0e-8;

    my $fac=2.0; my $sum=0; my $termbf=0; my $term=''; 

    my $a2 = -2.0*$alam*$alam;
    for (my $j=1; $j<=100; $j++) {
	$term = $fac*exp($a2*$j*$j);
	$sum += $term;
	return $sum if (abs($term) <= $EPS1*$termbf or abs($term)<=$EPS2*$sum);
	$fac = -$fac;
	$termbf = abs($term);
    }
    return 1.0;
}



sub randperm {
    use strict;
    use warnings;
    
    my ($N, $num) = @_;
    my @a = ();
    $#a = $N-1;
    for (my $i=1; $i<=$N; $i++) {
	$a[$i-1]=$i;
    }

    my $n = ''; my @output = (); $#output=$num-1;
    my $k=0;
    while ($num-->0) {
	$n = int(rand(scalar(@a)));
	
	$output[$k++] = $a[$n];
	splice @a, $n, 1;
    }

    return @output;
}


sub generate_rand_int {
    use strict;
    use warnings;

    my ($upper, $num) = @_;
    my @a = ();

    while ($num-->0) {
	push @a, int(rand($upper));
    }
    
    return @a;
}




#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
  
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    } else{
	return '*';
	#print STDERR "Bad codon \"$codon\"!!\n";
	#exit;
    }
}

sub seq2protein {
    my ($seq) = @_;
    
    my $protein = '';
    my $L = length($seq);
    for (my $j=0; $j<$L-2; $j+=3) {
	$aa = codon2aa(substr($seq, $j, 3));
	$protein .= $aa;
    }
    return $protein;
}




# search sorted_array for key
# return index of matching element if found
# otherwise return index of element if it is inserved into the array
sub binary_search {
    my ($key, @sortedArray) = @_;
    use strict;
    use warnings;
    
    my $first = 0;
    my $last = $#sortedArray;
    my $mid = ''; 

    while ($first <= $last) {
	$mid = sprintf('%d',($first + $last)/2);
	if ($key > $sortedArray[$mid]) {
	    $first = $mid + 1;
	} elsif ($key < $sortedArray[$mid]) {
	    $last = $mid - 1;
	} else {
	    return $mid;
	}
    }

    return $first;
}


# search sorted_array for key
# return index of matching element if found
# otherwise return index of element if it is inserved into the array
sub binary_search_ref {
    my ($key, $sortedArrayRef) = @_;
    use strict;
    use warnings;
    
    my $first = 0;
    my $last = $#{$sortedArrayRef};
    my $mid = ''; 

    while ($first <= $last) {
	$mid = sprintf('%d',($first + $last)/2);
	if ($key > ${$sortedArrayRef}[$mid]) {
	    $first = $mid + 1;
	} elsif ($key < ${$sortedArrayRef}[$mid]) {
	    $last = $mid - 1;
	} else {
	    return $mid;
	}
    }

    return $first;
}




sub pearson_xcorr {
    my ($refA, $refB) = @_;
    use strict;
    use warnings;
    
    my $N = scalar @{$refA};
    die "Error in $N != ", scalar(@{$refB}) if ($N ne scalar(@{$refB}));

    my ($mu, $sig) = ();
    $mu = mean(@{$refA});
    $sig = std(@{$refA});
    return 0 if ($sig <= 0);

    foreach (@{$refA}) {
	$_ = ($_-$mu)/$sig;
    }
    
    $mu = mean(@{$refB});
    $sig = std(@{$refB});
    return 0 if ($sig <= 0);

    foreach (@{$refB}) {
	$_ = ($_-$mu)/$sig;
    }

    my $d = 0; 
    for (my $i=0; $i<$N; $i++) {
	$d += $refA->[$i] * $refB->[$i];
    }

    return $d/$N;
}


sub spearman_rank_corr {
    my ($refA, $refB) = @_;
    use strict;
    use warnings;
    
    my $N = scalar @{$refA};
    die "Error in $N" if ($N ne scalar(@{$refB}));
    my $d = 0; my $tmp = 0;
    for (my $i=0; $i<$N; $i++) {
	$tmp = ($refA->[$i]-$refB->[$i]);
	$d += $tmp*$tmp;
    }
    return 1 - 6*$d/$N/($N*$N-1);
}


sub consensus2seqlist {
    my ($cons) = @_;
    use strict;
    use warnings;
    
    my @seqList = ();    
    my %iub2character_class = (
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
	N => '[ACGT]',
	n => '[ACGT]',		       
    );
    
    $cons = uc $cons;
    my $L = length($cons);  my @a=(); my $t=''; my $b='';  my @fooList=();
    for (my $i=0; $i<$L; $i++) {
	$t = $iub2character_class{substr($cons,$i,1)};
	if (length($t)==1) {
	    if (not @seqList) {
		push @seqList, $t;
	    } else {
		foreach (@seqList) {
		    $_ .= $t;
		}
	    }
	    next;
	}

	$t=~s/[\]\[]//g;
	@a = split("",$t);
	if (not @seqList) {
	    push @seqList, @a;
	    next;
	}

	@fooList=();
	foreach $b (@a) {
	    foreach (@seqList) {
		push @fooList, $_ . $b;
	    }
	}
	
	@seqList = @fooList;
    }
    return @seqList;
}




# check the number of matched bases between two sequences
# start all from position 0. 
# if caseSensitive=1, then without consider cases
sub check_num_match_between_two_seqs {
    my ($seqA, $seqB, $caseSensitive) = @_;
    use strict;
    use warnings;
    
    if ($caseSensitive != 1) {
	$seqA = uc $seqA; $seqB = uc $seqB;
    }
    
    my @listA = split("", $seqA);
    my @listB = split("", $seqB);
    
    my $L = min(scalar(@listA), scalar(@listB));
    my $n = 0;
    for (my $i=0; $i<$L; $i++) {
	$n++ if ($listA[$i] eq $listB[$i]);
    }
    return $n;
}


sub read_fa_to_seq_list {
    my ($file) = @_;
    use strict;
    use warnings;

    open myTmpIN, $file or die "cannot open $file";
    $/=">";
    $_ = <myTmpIN>;
    
    my @outList=(); my ($id,$seq)=();
    while (<myTmpIN>) {
	chomp;
	($id,$seq)=split (/\n/,$_,2);
	$seq =~ s/\n//g;
	push @outList, {'id'=>$id, 'seq'=>$seq};

    }
    close myTmpIN;
    $/="\n";
    return @outList;
}


# last line of perl module 
1;
