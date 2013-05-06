#!/usr/bin/env python
from optparse import OptionParser
import os
import sys
from Bio import SeqIO, motifs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from multiprocessing import Pool
from random import gauss
from numpy import array, cumsum, sort, sum, searchsorted, newaxis
from numpy.random import rand
from itertools import repeat, chain, izip

"""
Weighted random selection
returns n_picks random indexes.
the chance to pick the index i 
is give by the weight weights[i].
"""
def weighted_pick(weights,n_picks):
    t = cumsum(weights)
    s = sum(weights)
    return searchsorted(t,rand(n_picks)*s)
 
"""
A mappable function for generating a sequence. Generates
a DNA sequence according to the parameters.
Input:
pwm - the position weight matrix. 
bg - vector of background model
mu - the mean size of the sequence. >= rows of pwm
std - standard deviation of sequence size. if 0, uniform lengths

Output:
seq - the sequence string
"""
def genseq(bg, pwm, mu = 200, std = 0):
    if pwm == None:
        motifwidth = 0#if there is no motif, then 0 motif width
    else:
        motifwidth = pwm.shape[0]#width of the motif is number of rows in pwm
    if std == 0:
        seqlength = mu
    else:
        seqlength = max(motifwidth,int(gauss(mu,std)))#I don't want anything less than the motif width
    bgseqslength = seqlength - motifwidth
    bgseqind =  weighted_pick(bg,bgseqslength)#randomly select characters' indices
    letters = array(['A','C','G','T'])
    bgseq = list(letters[bgseqind])#background list
    if motifwidth > 0:
        randmotifseq = "".join([letters[weighted_pick(pwmrow,1)[0]] for pwmrow in pwm])
        bgseq.insert(bgseqslength/2,randmotifseq)
    seq = "".join(bgseq)#insert motif in middle and join string
    return seq

"""

Input:
fracsMotif - a list of fractions of the motif(s)
"""
def makeFASTA(fastafilename, bg, pwms, fracsMotif, n):
    output_handle = open(fastafilename, 'w')
    bgfrac = 1 - sum(fracsMotif)
    fracsMotif.insert(0,bgfrac)#first element is fraction of no motif
    fracsMotif = array(fracsMotif)#convert to array to allow multiplication
    freqsMotif = n*fracsMotif#number of each sequence type to generate
    freqsMotif = freqsMotif.round().astype(int)
    pwms.insert(0,None)#background acts as a the first pseudo pwm, so no motif first
    pwmiters = [repeat(pwm, freq) for pwm, freq in zip(pwms, freqsMotif)]#a list of iterables for sequence generation
    pwmchain = chain(*pwmiters)
    l = sum(freqsMotif)
    frags = [SeqRecord(Seq(genseq(bg,pwm),IUPAC.ambiguous_dna),'fragment_%i' % (i+1),'','') for pwm, i in izip(pwmchain,xrange(l))]
    #I just want to avoid for loops. Generate all seqrecords to be written to file
    SeqIO.write(frags, output_handle, "fasta")
    output_handle.close()

"""
Accepts the filename of a motif file, its type, and returns a list of numpy pwm arrays
"""
def motifParse(filename, filetype):
    if filename == None:
        return None#no filename, return empty matrix
    f = open(filename, 'r')
    s = f.read()
    f.close()
    
    if filetype == 'motiflist':#the format from xhx's motiflist
        t = s.split('\n')[1:-1]#split by line. first line is info. last line empty
        u = [a.split('\t')[:-1] for a in t]#the last character is consensus sequence
        pwm = array(u,dtype=float)
        pwm = pwm/pwm.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob
        pwms = [pwm]
    
    return pwms

"""
Accepts the filename of a background file, its type, and returns a numpy pwm array
"""
def backgroundParse(filename):
    if filename == None:
        return array([0.25,0.25,0.25,0.25])#default background distribution

"""
Accepts a string float and returns a float list. If float string cannot be converted,
it is assumed it is a filename. The filename should have the same number of lines as
there are motifs. 
"""
def fracsParse(fracstring):
    try:
        fracs = [float(fracstring)]
    except ValueError:
        fracsfile = open(fracstring,'r')
        r = fracsfile.read()
        rs = r.split('\n')[:-1]#split lines. ignore last line, which is empty
        fracs = [float(s) for s in rs]
        fracsfile.close()
    return fracs

if __name__ == "__main__":
    usage = "usage: %prog [options] <output FASTA>"
    description = "This program generates a FASTA file containing randomized sequences appropriate for analysis by OnlineMEME"
    parser = OptionParser(usage=usage,description=description)
    parser.add_option("-p", "--processes", help="optional number of parallelized processes")
    parser.add_option("-m", "--motif", help="File holding motif(s). Default: no motifs")
    parser.add_option("-b", "--background", help="Background file. Default: 0-order equal distribution")    
    parser.add_option("-n", "--numsequences", help="Number of sequences to write. Default:100", default="100")
    parser.add_option("-f", "--fraction", help="Fraction of sequences containing motif. Default:0.5", default="0.5")    
    (options, args) = parser.parse_args()
    #As of 4/9/13, work on motif parsers
    pwms = motifParse(options.motif,'motiflist')
    bg = backgroundParse(options.background)
    if len(args) == 0:
        parser.print_help()
        sys.exit()
    fastafilename = args[0]
    n = int(options.numsequences)
    fracsMotif = fracsParse(options.fraction)
    makeFASTA(fastafilename, bg, pwms, fracsMotif, n)
