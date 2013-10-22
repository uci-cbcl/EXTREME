README for EXTREME 1.0
========================

EXTREME is an efficient motif discovery algorithm. It applies the online EM algorithm to discover motifs. It uses the same
model as MEME, representing sequence binding preference as position frequency matricies (PFMs). EXTREME is written in Python,
and incorporates source code from MEME and DREME.


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
EXTREME comes in four different flavors:

* `EXTREME.py`. The original EXTREME algorithm. This flavor is best suited for ChIP-Seq datasets for TFs with long
motifs with high information content such as NRSF. Its seeding strategy involves generating starting points from
sequences centered around an instance of the DREME regular expression seed. It includes a preprocessing set
to take the middle 100 nucleotides of all sequences.

* `ShortEXTREME.py`. This flavor is best suited for ChIP-Seq datasets for TFs with 
short motifs such as NANOG. Its seeding strategy involves padding the DREME regular expression seed with universal
letters and then generating a PWM starting point from the nucleotide counts of sequences matching the regular expression.
It includes a preprocessing set to take the middle 100 nucleotides of all sequences.

* `FootprintEXTREME.py`. This flavor is best suited for finding motifs in footprint data. Uses the same seeding
strategy as ShortEXTREME, but pads using four universal letters on each side instead of two.

* `ParFootprintEXTREME.py`. This flavor is best suited for finding motifs in footprint data. It is meant to be
used as part of a batch job in parallel. It generates an initial PWM from a DREME regular expression seed in a
DREME output file and then runs the online EM algorithm to completion.

Arguments
---------
The following arguments are common to all variants of EXTREME:

* `-t TRIES`. The number of tries of the online EM algorithm before giving up on the current seed.
* `-m NUMMOTIFS`. The number of motifs to search for.
* `-o OUTPUT`. The output folder containing all output files.
* `-s SEED`. Random seed for shuffling sequences and dataset positions.
* `-w WIDTH`. The width of a motif to search for. If this value is set, both MINW and MAXW are set to this value.
* `-minw MINW`. The minimum motif width to search for.
* `-maxw MAXW`. The maximum motif width to search for. EXTREME.py looks for motifs of widths between MINW and MAXW,
in increments of sqrt(2). For the other three variants, if MAXW is greater than MAXK, it will attempt to generate
a new regular expression seed at width MAXW.

The following arguments are unique ShortEXTREME.py, FootprintEXTREME.py, and ParFootprintEXTREME.py:
* `-mink MINK`. The minimum DREME core width to search for.
* `-maxk MAXK`. The maximum DREME core width to search for. 

Running EXTREME
---------------
An example of running EXTREME using the included NRSF example:
```
$ python EXTREME.py -minw 15 -maxw 35 -m 5 -o outputFolder NRSF.fasta
```

Alternatively, the Python script file can be made an executable:
```
$ chmod +x EXTREME.py
$ ./EXTREME.py -minw 15 -maxw 35 -m 5 -o outputFolder NRSF.fasta
```

Output files
------------
**\*/Motif_x.png** PNG output of the x-th motif. Includes all motifs, not just the most significant ones.

**\*/Motif_x.eps** EPS output of the x-th motif. Includes all motifs, not just the most significant ones.

**\*/MEMEoutput.meme** Minimal MEME format output of discovered motifs.
