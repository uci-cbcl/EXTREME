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
from numpy import array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, power, ceil, prod
from numpy.random import rand
from itertools import repeat, chain, izip
from pygr import seqdb

"""
Equation 14 from the Bailey and Elkan paper. Calculates P(X|theta_motif). That is,
the probability of the sequence given the motif model. It is defined as the product
of the frequencies the letter at each position.

Input:
X, A string. The sequence. Must only contain A, C, G, T
theta_motif, PWM. Columns in order of A, C, G, T. Each row must sum to 1. Number of rows must be same as X length

Output:
p, A double. P(X|theta_motif) 
"""
def px_motif(X, theta_motif):
    d = {'A':0,'C':1,'G':2,'T':3}
    ps = [row[d[x]] for x, row in zip(X,theta_motif)]
    p = prod(ps)
    return p

"""
Equation 15 from the Bailey and Elkan paper. Calculates P(X|theta_background). That is,
the probability of the sequence given the background model. It is defined as the product
of the frequencies the letter at each position.

Input:
X, A string. The sequence. Must only contain A, C, G, T
theta_background, Array of probabilities for each nucleotide, in order of A, C, G, T. Must sum to 1

Output:
p, A double. P(X|theta_background) 
"""
def px_background(X, theta_background):
    d = {'A':0,'C':1,'G':2,'T':3}
    indices = [d[x] for x in X]
    ps = theta_background[indices]
    p = ps.prod()
    return p

"""
Equation 10 from the Bailey and Elkan paper. The expectation of the Z value for the sequence X,
given the motif frequency and the motif and background models. Note that g=2, where group 1
is the motif, and group 2 is the background.

Input:
X, a string. The sequence.
theta_motif, a matrix. The pwm of the motif.
theta_background, an array. The background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Notes:
5/6/13, work on this part
"""
def Z0(X,theta_motif, theta_background,lambda_motif):
    return px_motif(X,theta_motif)*lambda_motif/ \
        (px_motif(X,theta_motif)*lambda_motif + px_background(X,theta_background)*(1-lambda_motif))

def EM():
    #E-step
    #update the z's
    #M-step
    return 0

"""
Function that returns a list of subsequences, X,
of width W, from the dataset of sequences, Y.
It is effectively a ragged array

Input:
Y, dataset of sequences
W, motif width

Notes:
For the Online version, this may need to be removed for RAM purposes.
"""
def getSubsequences(Y, W):
    return [[y[k:k+W] for k in xrange(len(y)-W+1)] for y in Y]

"""
The main MEME algorithm.

Input:
Y, dataset of sequences
W, width of motifs to search for
NPASSES, number of distinct motifs to search for

5/14/13:
Assume that each sequence is the size of the motif for now.
That is, Y = X.
"""
def meme(Y,W,NPASSES):
    X = getSubsequences(Y)#this step may need to be moved for Online to save RAM
    n = sum([len(x) for x in X])#number of subsequences. this step may need to be changed for Online
    N = len(Y)#number of sequences
    alpha = 0.99#probability of catching sequence with motif
    for i in range(NPASSES):#search for NPASSES of motifs of width W
        a = arange(0,log2(n/(2*W*sqrt(N)))+1,1)#the +1 ensures the last guess included
        b = a + log2(sqrt(N)/n) 
        lambda0s = power(2,b)#array holding the heuristic initial lambda values to try
        for lambda0 in lambda0s:
            Q = log2(1-alpha)/log2(1-lambda0)#number of sequences to sample to catch motif with prob alpha
            for j in sample(xrange(n),ceil(Q)):#guarantee to sample at least one sequence with motif
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