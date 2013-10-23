README for EXTREME 1.0
========================

EXTREME is an efficient motif discovery algorithm. It applies the online EM algorithm to discover motifs. It uses the same
model as MEME, representing sequence binding preference as position frequency matricies (PFMs). EXTREME is written in Python,
and incorporates source code from MEME and DREME (by T. Bailey) and are part of the [MEME Suite](http://meme.nbcr.net/meme/)


INSTALL
=======

Prerequisites
-------------
* Python (2.7). [Python 2.7.3](http://www.python.org/download/releases/2.7.3/) is recommended.

* [Numpy](http://www.numpy.org/)(>=1.6.1).

* [WebLogo] (https://code.google.com/p/weblogo/) (3).

* [Cython] (http://cython.org).

* [gcc] (http://gcc.gnu.org/).

Although not required by EXTREME, [bedtools] (https://github.com/arq5x/bedtools/) can be useful for manipulating BED 
files and generating FASTA files.

Install from source
-------------------
Download the compressed source files, cd into the folder, and use the following command to compile the Cython bindings 
to the MEME source files:

```
$ python setup.py build_ext --inplace
```


USAGE
=====


Overview
--------
EXTREME comes in three different flavors:

* `EXTREME.py`. The original EXTREME algorithm. This flavor is best suited for ChIP-Seq datasets for TFs with long
motifs with high information content such as NRSF. It uses the first seeding strategy which involves generating starting points from
sequences centered around an instance of the DREME regular expression seed. It includes a preprocessing step
to trim all sequences to the middle 100 bps.

* `ShortEXTREME.py`. This flavor is best suited for ChIP-Seq datasets for TFs with 
short motifs such as NANOG or CTCF. Its seeding strategy involves padding the DREME regular expression seed with universal
letters and then generating a PFM starting point from the nucleotide counts of sequences matching the regular expression.
It includes a preprocessing set to take the middle 100 nucleotides of all sequences.

* `ParallelFootprintEXTREME.py`. This flavor is best suited for finding motifs in footprint data. It uses the second seeding strategy.
PFMs are initialized from JASPAR or ENCODE instead of generating one using DREME. We have supplied JASPAR and ENCODE
databases in Minimal MEME Format for the seeding. It is parallelizable, processing each seed independently with the
online EM algorithm.


Arguments
---------
The following arguments are common to all variants of EXTREME:

* `-t TRIES`. The number of tries of the online EM algorithm before giving up on the current seed.
* `-s SEED`. Random seed for shuffling sequences and dataset positions.

The following arguments are common to EXTREME.py and ShortEXTREME.py:
* `-m NUMMOTIFS`. The number of motifs to search for.
* `-o OUTPUT`. The output folder containing all output files.
* `-w WIDTH`. The width of a motif to search for. If this value is set, both MINW and MAXW are set to this value.
* `-minw MINW`. The minimum motif width to search for.
* `-maxw MAXW`. The maximum motif width to search for. EXTREME.py looks for motifs of widths between MINW and MAXW,
in increments of sqrt(2). For ShortEXTREME.py, if MAXW is greater than MAXK, it will attempt to generate
a new regular expression seed at width MAXW using an experimental long regular expression generator in DREME.

The following arguments are unique ShortEXTREME.py:
* `-mink MINK`. The minimum DREME core width to search for.
* `-maxk MAXK`. The maximum DREME core width to search for. 
* `-p PADDING`. The number of universal letters to add to each side of the regular expression seed.

The following arguments are unique ParallelEXTREME.py:
* `-p PSEUDOCOUNTS`. Uniform pseudo counts to add to initial PFM guess (default 0.1).

Running EXTREME
---------------
An example of running EXTREME using the included NRSF example. cd into the ExampleFiles directory and run the EXTREME
algorithm. This particular example searches for 2 motifs of width 23 and width 33 using 20 seeds per width. Results are put 
into the folder "outputFolder":
```
$ cd ExampleFiles
$ python ../EXTREME.py -minw 23 -maxw 33 -m 2 -o outputFolder -t 20 GM12878_NRSF_intersected.fasta.masked
```
The ../ indicates to execute EXTREME.py in the parent folder.

Alternatively, the Python script file can be made an executable:
```
$ chmod +x EXTREME.py
$ cd ExampleFiles
$ python ../EXTREME.py -minw 23 -maxw 33 -m 2 -o outputFolder -t 20 GM12878_NRSF_intersected.fasta.masked
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
> ../ParallelFootprintEXTREME.py -p 0.0 K562_footprints_extended5_merged_Round2.fasta.masked K562_footprints_extended5_merged_Round2_negative.fasta.masked ENCODE_Motif_names ENCODE_PosSites_Round2 ENCODE.meme $i
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
> ../ParallelFootprintEXTREME.py -p 0.0 K562_footprints_extended5_merged.fasta.masked K562_footprints_extended5_merged_negative.fasta.masked ENCODE_Motif_names ENCODE_PosSites ENCODE.meme $i
> EOF
> done
```

Note that these last two examples will generate many output folders. You have been warned.

Output files
------------
**\*/Motif_x.png** PNG output of the x-th motif. Includes all motifs, not just the most significant ones (that is, the final
result after convergence of any seed).

**\*/Motif_x.eps** Same as above, except in EPS format.

**\*/MEMEoutput.meme** Minimal MEME format output of discovered motifs. 
