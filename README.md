README for EXTREME 2.0
========================

EXTREME is an efficient motif discovery algorithm. It applies the online EM algorithm to discover motifs. It uses the same
model as MEME, representing sequence binding preference as position frequency matricies (PFMs). EXTREME is written in Python,
and incorporates source code from MEME and DREME (by T. Bailey), which are part of the [MEME Suite](http://meme.nbcr.net/meme/).
A paper for EXTREME has been submitted to the journal Bioinformatics: Oxford Journals and is currently under peer review.


INSTALL
=======

Required
--------
* [Python 2.7](http://www.python.org/download/releases/2.7.6/).

* [Numpy](http://www.numpy.org/).

* [Perl] (http://www.perl.org/).

* [Java] (http://www.java.com/).

Optional
________

* [WebLogo] (https://code.google.com/p/weblogo/) (3). For generating visual sequence logos in eps and png format.

* [Cython] (http://cython.org). For building C bindings to the MEME source code. This is required if you want to calculate E-values.

* [gcc] (http://gcc.gnu.org/). For compiling C code. Cython needs this.

* [bedtools] (https://github.com/arq5x/bedtools/). Can be useful for manipulating BED files and generating FASTA files. When generating your own FASTA files, we recommend using the masked reference genomes to extract genomic sequences from.

* [MEME suite] (http://meme.nbcr.net/meme/). Has a bunch of useful tools for motifs.


Install from source
-------------------
Download the latest release (zip (https://github.com/uci-cbcl/EXTREME/archive/2.0.0.zip) tar.gz (https://github.com/uci-cbcl/EXTREME/archive/2.0.0.tar.gz)) and decompress. 

Optional: If you want to calculate E-values for your motifs, then you need to build Cython bindings to the MEME source files. Keep in mind that Cython and gcc are usually difficult to work with. I have had the best success on a Linux setup. cd into the src folder, and use the following command:

```
$ python setup.py build_ext --inplace
```


USAGE
=====

Arguments
---------

The following are arguments for GappedKmerSearch.py, the word searching algorithm for the seeding:
* `-t TRIES`. The number of different bias factors to try before giving up on the current seed.
* `-s SEED`. Random seed for shuffling sequences and dataset positions.
* `-p PSEUDOCOUNTS`. Uniform pseudo counts to add to initial PFM guess (default 0.0).
* `-minsites MINSITES`. Minimum number of sites the motif should have (default 10).

The following are arguments 


The following are arguments for EXTREME.py, the EXTREME algorithm:

* `-t TRIES`. The number of different bias factors to try before giving up on the current seed.
* `-s SEED`. Random seed for shuffling sequences and dataset positions.
* `-p PSEUDOCOUNTS`. Uniform pseudo counts to add to initial PFM guess (default 0.0).
* `-minsites MINSITES`. Minimum number of sites the motif should have (default 10).
* `-maxsites MAXSITES`. Minimum number of sites the motif should have. If not specified, it is set to five times the number of predicted motif sites based on the initial PFM guess
* `-saveseqs SAVESEQS`. A switch. If used, the positive and negative sequence set will be saved to Positive_seq.fa and Negative_seq.fa, respectively, with instances of the discovered motif replaced with capital Ns.

Running EXTREME
---------------
An example of running EXTREME using the included NRSF example. cd into the ExampleFiles directory. First, we need to generate some seeds:
```
$ python ../GappedKmerSearch.py -l 8 -ming 33 -m 2 -o outputFolder -t 20 GM12878_NRSF_intersected.fasta.masked
```


```
$ cd ExampleFiles
$ python ../EXTREME.py -minw 23 -maxw 33 -m 2 -o outputFolder -t 20 GM12878_NRSF_intersected.fasta.masked
```
The ../ indicates to execute EXTREME.py in the parent folder. Note that EXTREME is fast enough to test many seeds
and process them to convergence with the online EM algorithm. However, a more thorough search using more seeds
will linearly increase the running time of EXTREME.

Alternatively, the Python script file can be made an executable:
```
$ chmod +x EXTREME.py
$ cd ExampleFiles
$ ../EXTREME.py -minw 23 -maxw 33 -m 2 -o outputFolder -t 20 GM12878_NRSF_intersected.fasta.masked
```
However, this approach is not always successful. Most motifs are usually shorter and do not work so well with this strategy, such
as CTCF and NANOG. In cases like these, we recommend using ShortEXTREME.py. Here is an example that uses ShortEXTREME.py
to discover a motif in the included CTCF dataset. It uses DREME's experimental long regular expression algorithm to discover
a motif of width 12, and then pads the motif on both sides with 2 universal letters to generate an initial seed of 
width 16. Results are put into the folder "outputFolder":
```
$ cd ExampleFiles
$ python ../ShortEXTREME.py -minw 3 -maxw 12 -p 2 -m 1 -o outputFolder SKNSHRA_CTCF_intersected.fasta.masked
```
Now, this is probably what you came here for. How can you discover a bunch of motifs in footprints simultaneously
very quickly? First, let's start with searching for one motif. Here is an example where we use the included JASPAR
database to seed with the CTCF motif.
```
$ cd ExampleFiles
$ ../ParallelFootprintEXTREME.py -p 0.1 K562_footprints_extended5_merged.fasta.masked K562_footprints_extended5_merged_negative.fasta.masked JASPAR_Motif_names JASPAR_PosSites JASPAR_CORE_2009.meme 156
```
Unlike the other two flavors, this last flavor has 5 required arguments. The first one is the dataset in FASTA format.
The second one is a "negative dataset", which we generated by shuffling dinucleotides by using the sequence.py script. We
require a negative dataset so that users can compare motif instances in positive and negative datasets to decide whether
a motif is valid.
The third one is a file containing the motif names in order. The fourth one is a file that lists the number of sites
discovered using the initial guess. Users can generate their own file like this by using FIMO with a Minimal MEME file
and a FASTA file. The last argument is an integer, which tells the script which motif to use. The first FASTA file was 
generated by taking K562 footprints from the ENCODE footprint dataset, extending each footprint by 5 bp, and merging
all intersecting regions. Unlike the other two variants, you do not get a choice where results are outputted to. Results
are always put into a folder that shares the name of the motif seed used.

Now, if you want to parallelize this, you will need a cluster with a lot of cores and memory. Here is an example
in which we go through all 683 seeds in the ENCODE database:
```
$ cd ExampleFiles
$ for i in {1..683}; do
> qsub -S /bin/bash -q mycluster <<EOF
> ../ParallelFootprintEXTREME.py -t 15 -p 0.0 K562_footprints_extended5_merged_Round2.fasta.masked K562_footprints_extended5_merged_Round2_negative.fasta.masked ENCODE_Motif_names ENCODE_PosSites_Round2 ENCODE.meme $i
> EOF
> done
```
Note that in this example, we used 0 pseudo counts, which differs from the default. The "Round2" refers to how we took
the FASTA file from the previous example and deleted four motif instances and repetitive elements so that the
EXTREME algorithm can converge onto other motifs.

Here's another example for your amusement:
```
$ cd ExampleFiles
$ for i in {1..683}; do
> qsub -S /bin/bash -q mycluster <<EOF
> ../ParallelFootprintEXTREME.py -t 15 -p 0.0 K562_footprints_extended5_merged.fasta.masked K562_footprints_extended5_merged_negative.fasta.masked ENCODE_Motif_names ENCODE_PosSites ENCODE.meme $i
> EOF
> done
```

Note that these last two examples will generate many output folders. You have been warned.

Output files
------------
**\*/Motif_x.png** PNG output of the x-th motif. Includes all motifs, not just the most significant ones (that is, the final
result after convergence of any seed).

**\*/Motif_x.eps** Same as above, except in EPS format.

**\*/MEMEoutput.meme** Minimal MEME format output of discovered motifs (not all seeds. Only the motifs EXTREME selected at the end
of a seed search.)
