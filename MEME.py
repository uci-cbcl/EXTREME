#!/usr/bin/env python
from optparse import OptionParser
import os
import sys
from Bio import SeqIO, motifs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from multiprocessing import Pool
from random import gauss, sample
from numpy import array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, power, ceil
from numpy.random import rand
from itertools import repeat, chain, izip
from pygr import seqdb

def EM():
    return 0

def meme():
    Y = 0#dataset of sequences
    W = 0#width of motifs to search for
    NPASSES = 0#number of distinct motifs to search for
    X = 0
    n = len(X)#number of subsequences
    N = len(Y)#number of sequences
    alpha = 0.99#probability of catching sequence with motif
    for i in range(NPASSES):#search for NPASSES of motifs of width W
        a = arange(0,log2(n/(2*W*sqrt(N))),1)
        b = a + log2(sqrt(N)/n) 
        lambda0s = power(2,b)#array holding the heuristic initial lambda values to try
        for lambda0 in lambda0s:
            Q = log2(1-alpha)/log2(1-lambda0)#number of sequences to sample to catch motif with prob alpha
            for j in sample(xrange(n),ceil(Q)):#guarantee to sample at least one sequence
                seq = X[j]#pick subsequence j
    

if __name__ == "__main__":
    usage = "usage: %prog [options] <input FASTA>"
    description = "This program generates a FASTA file containing randomized sequences appropriate for analysis by OnlineMEME"
    parser = OptionParser(usage=usage,description=description)
    parser.add_option("-p", "--processes", help="optional number of parallelized processes")
    parser.add_option("-m", "--motif", help="File holding motif(s). Default: no motifs")
    parser.add_option("-b", "--background", help="Background file. Default: 0-order equal distribution")    
    parser.add_option("-n", "--numsequences", help="Number of sequences to write. Default:100", default="100")
    parser.add_option("-f", "--fraction", help="Fraction of sequences containing motif. Default:0.5", default="0.5")    
    (options, args) = parser.parse_args()