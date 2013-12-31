#!/usr/bin/perl -w
use FindBin qw($Bin);
use lib "$Bin";
use gene_utils;


die "usage: perl $0 consensus_file" if ($#ARGV<0);






%code2prob = ();
push @{$code2prob{'A'}},  (1, 0, 0, 0);
push @{$code2prob{'C'}}, (0, 1, 0, 0);
push @{$code2prob{'G'}}, (0, 0, 1, 0);
push @{$code2prob{'T'}}, (0, 0, 0, 1);
push @{$code2prob{'R'}}, (0.5, 0, 0.5, 0);
push @{$code2prob{'Y'}}, (0, 0.5, 0, 0.5);
push @{$code2prob{'K'}}, (0, 0, 0.5, 0.5);
push @{$code2prob{'M'}}, (0.5, 0.5, 0, 0);
push @{$code2prob{'S'}}, (0, 0.5, 0.5, 0);
push @{$code2prob{'W'}}, (0.5, 0, 0, 0.5);
push @{$code2prob{'N'}}, (0.25, 0.25, 0.25, 0.25);



@all = get_file_data($ARGV[0]); chomp @all;




foreach (@all) {
    chomp $_;
    @a = split(" ");
    $cons = $a[0];
    $cons = uc $cons;
    print ">$cons\n";
    
    @a = split("", $cons);
    foreach $b (@a) {
	@prob = @{$code2prob{$b}};
	foreach (@prob) {
	    $_ = ($_+0.01)/1.04;
	}
	$s = sum(@prob);
	foreach (@prob) {
	    $_ = sprintf('%5.4f',$_/$s);
	}

	print join("\t",@prob),"\t$b", "\n";
    }

    print "\n\n";

	
    
    

}
