README for EXTREME 2.0
========================

NOTE: As of 2018, this repository is deprecated. Unless you have specific reason to use EXTREME, I suggest to use [YAMDA](https://github.com/daquang/YAMDA),

EXTREME is an efficient motif discovery algorithm. It applies the online EM algorithm to discover motifs. It uses the same
model as MEME, representing sequence binding preference as position frequency matricies (PFMs). EXTREME is written in Python,
and incorporates source code from MEME and DREME (by T. Bailey), which are part of the [MEME Suite](http://meme.nbcr.net/meme/).
EXTREME is published. You can download an Advance Access version of the paper [here](https://github.com/uci-cbcl/EXTREME/blob/master/AdvanceAccess.pdf). If you have any more questions you can e-mail me at dxquang@uci.edu.


Citing EXTREME
========================

Quang, D., & Xie, X. (2014). EXTREME: an online EM algorithm for motif discovery. Bioinformatics, btu093.

INSTALL
=======

Required
---------
* [Python >= 2.7](http://www.python.org/download/releases/2.7.6/). The majority of the EXTREME algorithm is written in Python. Although these libraries are standard with Python >= 2.7 releases, make sure you have the following libraries available on your setup: argparse, collections, random, copy, errno, sys, itertools, string, math, re, and time. 

* [Numpy](http://www.numpy.org/). A standard external Python package for scientific computing.

* [Perl] (http://www.perl.org/).

* [Java] (http://www.java.com/).

Optional
--------

* [WebLogo] (https://weblogo.googlecode.com/files/weblogo-3.3.tar.gz) (3.3). For generating visual sequence logos in eps and png format. Whatever you do, do not download 3.4!

* [Cython] (http://cython.org). For building C bindings to the MEME source code. This is required if you want to calculate E-values.

* [gcc] (http://gcc.gnu.org/). For compiling C code. Cython needs this.

* [bedtools] (https://github.com/arq5x/bedtools/). Can be useful for manipulating BED files and generating FASTA files. When generating your own FASTA files, we recommend using the masked reference genomes to extract genomic sequences from.

* [MEME suite] (http://meme.nbcr.net/meme/). Has a bunch of useful tools for motifs.


Install from source
-------------------
Download the latest release ([zip] (https://github.com/uci-cbcl/EXTREME/archive/v2.0.0.zip) [tar.gz] (https://github.com/uci-cbcl/EXTREME/archive/v2.0.0.tar.gz)) and decompress. 

Optional: If you want to calculate E-values for your motifs, then you need to build Cython bindings to the MEME source files. Keep in mind that Cython and gcc are usually difficult to work with. I have had the best success on a Linux setup. cd into the src folder, and use the following command:

```
$ python setup.py build_ext --inplace
```


USAGE
=====

Arguments
---------

The following are arguments for GappedKmerSearch.py, the word searching algorithm for the seeding:
* `-l HALFLENGTH`. The number of exact letters on each side of the word (default 4).
* `-ming MINGAP`. The minimum number of universal wildcard letters in the middle (default 0).
* `-maxg MAXGAP`. The maximum number of universal wildcard letters in the middle (default 10).
* `-minsites MINSITES`. Minimum number of sites a word should have to be included (default 10).
* `-zthresh ZTHRESHOLD`. Minimum normalized z-score for a word to be saved. A lower threshold increases the number of words saved (default 5).

The following are arguments for run_consensus_clusering_using_wm.pl, the hierarchical clustering algorithm for the seeding:
* `THRESHOLD`. The threshold for the clustering. Has values between 0 and 1. A value closer to 1 decreases the number of clusters, while a value closer to 0 increases the number of clusters. Recommended value is 0.3.


The following are arguments for EXTREME.py, the EXTREME algorithm:

* `-t TRIES`. The number of different bias factors to try before giving up on the current seed.
* `-s SEED`. Random seed for shuffling sequences and dataset positions.
* `-p PSEUDOCOUNTS`. Uniform pseudo counts to add to initial PFM guess (default 0.0).
* `-q INITALSTEP`.  The initial step size for the online EM algorithm. A VERY sensitive parameter. I get best success for ChIP size data (about 100,000 to 1,000,000 bps) with a step size of 0.05. For DNase footprinting, which usually has >5,000,000 bps, I find 0.02 works best (default 0.05).
* `-minsites MINSITES`. Minimum number of sites the motif should have (default 10).
* `-maxsites MAXSITES`. Minimum number of sites the motif should have. If not specified, it is set to five times the number of predicted motif sites based on the initial PFM guess
* `-saveseqs SAVESEQS`. A switch. If used, the positive and negative sequence set will be saved to Positive_seq.fa and Negative_seq.fa, respectively, with instances of the discovered motif replaced with capital Ns.
* `-b BACKGROUND`. A switch. If used, the minimal MEME output will use the background frequencies from the learning, instead of the default uniform frequencies.

Running EXTREME
---------------
An example of running EXTREME using the included ENCODE GM12878 NRSF ChIP-Seq dataset. cd into the ExampleFiles directory. First, we need to generate some seeds:
```
$ python ../src/fasta-dinucleotide-shuffle.py -f GM12878_NRSF_ChIP.fasta > GM12878_NRSF_ChIP_shuffled.fasta
$ python ../src/GappedKmerSearch.py -l 8 -ming 0 -maxg 10 -minsites 5 GM12878_NRSF_ChIP.fasta GM12878_NRSF_ChIP_shuffled.fasta GM12878_NRSF_ChIP.words
$ perl ../src/run_consensus_clusering_using_wm.pl GM12878_NRSF_ChIP.words 0.3
$ python ../src/Consensus2PWM.py GM12878_NRSF_ChIP.words.cluster.aln GM12878_NRSF_ChIP.wm
```
The first line generates a dinucleotide shuffled version of the positive sequence set to serve as a negative sequence set. The second line finds gapped words with two half-sites of length 8, between 0 and 10 universal wildcard gap letters, and at least 5 occurrences in the positive sequence set. The third line clusters the words and outputs the results to GM12878_NRSF_ChIP.words.cluster.aln (run_consensus_clusering_using_wm.pl always outputs results to the input filename with ‘cluster.aln’ appended at the end). The last line converts the clusters into PFMs which can be used as seeds for the online EM algorithm. These PFMs are saved in GM12878_NRSF_ChIP.wm. For your own data, you may need to play around with the parameters to get a good set of seeds.

Now let’s run the online EM algorithm.
```
$ python ../src/EXTREME.py GM12878_NRSF_ChIP.fasta GM12878_NRSF_ChIP_shuffled.fasta GM12878_NRSF_ChIP.wm 1
```

EXTREME.py uses PFM seeds from GM12878_NRSF_ChIP.wm to initialize the online EM algorithm. The last argument tells EXTREME which of these seeds to use. GM12878_NRSF_ChIP.wm should have 23 PFM seeds, so the last argument can be any value between 1 and 23 in this case. 

We have also included an ENCODE K562 DNase-Seq dataset. Try running EXTREME on your own with this dataset. In our publication, we used the parameters l=4, ming=0, maxg=10, minsites=10, zthresh=5 for the word search portion of the seeding. We also used an initial step size of q=0.02. You can imagine the initial step size as a sort of "shaking" parameter. A larger initial step corresponds to a more vigorous shaking, while a smaller value corresponds to a more gentle shaking. You can try experimenting with other sets of parameters too. Please keep me updated on what you find.

Output files
------------
EXTREME.py outputs files to a folder with the same name as seed the online EM algorithm is initialized from. For example, the first seed in
our NRSF example has the name “cluster1”, so all files will be output to the “cluster1” folder.

**\*/Motif_x.png** PNG output of the x-th motif. Includes all motifs, not just the most significant ones (that is, the final
result after convergence of any seed).

**\*/Motif_x.eps** Same as above, except in EPS format.

**\*/MEMEoutput.meme** Minimal MEME format output of discovered motifs (not all seeds. Only the motifs EXTREME selected at the end
of a seed search.)
