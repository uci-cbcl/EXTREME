#!/usr/bin/perl -w
use FindBin qw($Bin);
use lib "$Bin";
use gene_utils;


die "usage: perl name foo.motif.list foo.cluster.index" if ($#ARGV<1);

@motifList = get_file_data($ARGV[0]); chomp @motifList;

@outList = ();

open IN, $ARGV[1] or die "cannot open $ARGV[1]";
while (<IN>) {
    next if (not $_=~/^>/);
    @a = ();
    while (<IN>) {
	last if (not $_=~/[a-z0-9]/i);
	chomp;
	push @a, $motifList[$_];
    }
    push @outList, {'num'=>scalar(@a), 'list'=>join("\n",@a)};
}
close IN;
    

@outList = sort {$b->{'num'} <=> $a->{'num'}} @outList;

$k = 0;
foreach (@outList) {
    $k++;
    print ">cluster$k\t", $_->{'num'},"\n";
    print $_->{'list'},"\n\n";
}


