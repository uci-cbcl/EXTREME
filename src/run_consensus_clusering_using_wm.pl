#!/usr/bin/perl -w
use FindBin qw($Bin);
use lib "$Bin";
use gene_utils;

die "usage: perl name foo.consensus.list th_score([0,1])" if ($#ARGV<1);

$f = $ARGV[0];
$th_score = $ARGV[1];
$missing_penalty = 1;

$currDir = $Bin;
 
system "perl $currDir/consensus2weightmatrix.pl $f > $f.wm";

system "java -Xmx300M -cp $currDir/motif.jar motif.GenerateSimMatrixForNewFormatWM $f.wm 0 $f.wm.sim";

@all = get_file_data("$f.wm.sim");
@outList=();
foreach (@all) {
    @a = split(" ");
    $a[5] = ($a[5]>0 ? '+':'-');
    $a[4] = 1-$a[4];
    push @outList, join("\t", @a[0,2,4,5,6]);
}
write_list_to_file("$f.wm.dist", @outList);

# using java based hierahcical clustering
system "java -Xmx2500m -cp $currDir/motif.jar motif.HierarchicalClustering $f.wm.dist $th_score $missing_penalty > $f.cluster";

# use JAVA clustering program, also based on group average
system "perl $currDir/update_java_out_cluster.pl $f $f.cluster > $f.cluster.st";
system "perl $currDir/align_sequence_in_cluster.pl $f.wm.dist $f $f.cluster.st > $f.cluster.aln";


unlink "$f.wm";
unlink "$f.wm.sim";
unlink "$f.wm.dist";
unlink "$f.cluster";
unlink "$f.cluster.st";
