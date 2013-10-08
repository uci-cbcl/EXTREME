README for EXTREME 1.0
========================

An online implementation of the MEME algorithm


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

Output files
------------
