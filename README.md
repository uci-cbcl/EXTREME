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
algorithm. This particular example searches for 2 motifs of width 23 and width 33 using 20 seeds per width:
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

Output files
------------
**\*/Motif_x.png** PNG output of the x-th motif. Includes all motifs, not just the most significant ones.

**\*/Motif_x.eps** EPS output of the x-th motif. Includes all motifs, not just the most significant ones.

**\*/MEMEoutput.meme** Minimal MEME format output of discovered motifs.
