#!/usr/bin/env python
from optparse import OptionParser
import os
import sys
from Bio import SeqIO, motifs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from multiprocessing import Pool
from random import gauss, sample, shuffle
from numpy import mean,load,save,inf, sign, dot, diag, array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, log, power, ceil, prod, zeros, ones, concatenate
from numpy.random import rand
from pylab import imread, imshow, plot, show
from weblogolib import LogoData, LogoOptions, LogoFormat, png_formatter, unambiguous_dna_alphabet
from itertools import repeat, chain, izip
from pygr import seqdb
from bisect import bisect_left

"""
Equation 14 from the Bailey and Elkan paper. Calculates P(X|theta_motif). That is,
the probability of the sequence given the motif model. It is defined as the product
of the frequencies of the letter at each position.

This is a modified version that uses the indicator matrix instead. Saves the trouble
of calculating the indicator matrix at each step.

Input:
I, A boolean matrix. A representation of the sequence. Has dimension Wx4, the same as the PWM
theta_motif, PWM. Columns in order of A, C, G, T. Each row must sum to 1. Number of rows must be same as X length

Output:
p, A double. P(X|theta_motif) 
"""
def pI_motif(I, theta_motif):
    ps = theta_motif[I]#some fancy indexing tricks. Gets me an array of the relevant frequencies
    p = ps.prod()
    return p

"""
Equation 15 from the Bailey and Elkan paper. Calculates P(X|theta_background). That is,
the probability of the sequence given the background model. It is defined as the product
of the frequencies the letter at each position.

This is a modified version that uses the indicator matrix instead. Saves the trouble
of calculating the indicator matrix at each step. Note also that theta_background
has been switched to theta_background_matrix, essentially a PWM for the background.

Input:
I, A boolean matrix. A representation of the sequence. Has dimension Wx4, the same as the PWM
theta_background_matrix, Matrix of probabilities for each nucleotide, in order of A, C, G, T. 
Each row must sum to 1. This is basically like the usual array, but the row vector is repeated W times

Output:
p, A double. P(X|theta_background) 
"""
def pI_background(I, theta_background_matrix):
    ps = theta_background_matrix[I]#some fancy indexing tricks. Gets me an array of the relevant frequencies
    p = ps.prod()
    return p


"""
Equation 10 from the Bailey and Elkan paper. The expectation of the Z value for the sequence X,
given the motif frequency and the motif and background models. Note that g=2, where group 1
is the motif, and group 2 is the background. We can get the expectation of the Z value of 
the background by taking 1 - Z0.


This is a modified version that uses the indicator matrix instead. Saves the trouble
of calculating the indicator matrix at each step. 

Input:
I, indicator matrix. Represents a sequence.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
Z0 - Expected value of Z for the the indicator matrix I
"""
def Z0_I(I,theta_motif, theta_background_matrix,lambda_motif):
    a = pI_motif(I,theta_motif)*lambda_motif#saves a calculation
    b = pI_background(I,theta_background_matrix)*(1-lambda_motif)#saves another calculation
    Z0 = a/(a + b)
    return Z0

"""
Function that accepts a string, X, and returns an indicator matrix, I. Representing
each sequence as an indicator matrix will save time and memory.

Input:
X, a string. The sequence to be converted to an indicator matrix

Output:
I, a boolean matrix. The indicator matrix. Has dimensions Wx4. Each row is a position along the string.
Each column represents a nucleotide, in alphabetical order.
"""
def sequenceToI(X):
    l = array(list(X))#convert the string to an array of characters
    d = array(['A','C','G','T'])#the 4 nucleotides
    I = l[:,newaxis] == d#construct the matrix
    return I

"""
Computes a term in the summation, equation 11.

This is a modified version that uses the indicator matrix instead. Saves the trouble
of calculating the indicator matrix at each step. 

Input:
I, indicator matrix. Represents a sequence.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
elogL - a summation term of equation 11
"""
def expLog_I(I,theta_motif, theta_background_matrix,lambda_motif):
    a = pI_motif(I,theta_motif)*lambda_motif#saves a calculation
    b = pI_background(I,theta_background_matrix)*(1-lambda_motif)#saves another calculation
    Z0 = a/(a + b)
    elogL = Z0*log(a) + (1-Z0)*log(b)
    return elogL

"""
Calculates the expected log likelihood. Accepts the list of indicator matrices I, the current motif
frequency lambda_motif, and both PWMs. Returns the expected log likelihood.

Note that this uses the indicator matrices only in order to save time converting.
This function is only meant for debugging purposes, since in practice calculating the expected log likelihood
defeats the purpose of an online algorithm.

Input:
I, a list of indicator matrices. This might get confusing with I being a list here and a matrix elsewhere.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
expected_LogLikelihood, the expected log likelihood

7/5/13: Work on this next
"""
def expected_LogLikelihood(I, theta_motif, theta_background_matrix, lambda_motif):
    expected_LogLikelihood = 0
    for Ii in I:
        for Iij in Ii:
            expected_LogLikelihood = expected_LogLikelihood + expLog_I(Iij,theta_motif, theta_background_matrix,lambda_motif)
    return expected_LogLikelihood

"""
Absolute Euclidean distance between two arrays, u and v. This function
is for the EM algorithm's convergence.

Input:
u and v, arrays.

Output:
Euclidean distance between u and v.
"""
def dist(u,v):
    w = u - v
    w = power(w,2)
    return sqrt(w.sum())

"""
Function that returns a list of subsequences, X,
of width W, from the dataset of sequences, Y.
It is effectively a ragged array

Input:
Y, dataset of sequences. A pygr SequenceFileDB object 
W, motif width

Notes:
For the Online version, this may need to be removed for RAM purposes.
"""
def getSubsequences(Y, W):
    return [[str(Y[y][k:k+W]) for k in xrange(len(Y[y])-W+1)] for y in Y]

"""
A function to account for repetitive elements. Counts the number
of overlapping W-mers at each site on the sequence.

Input:
X, a string. The nucleotide sequence
W, motif width

Output:
counts, array of overlap occurrences. Length of L - W + 1
"""
def overlap_occurrences(X, W):
    L = len(X)
    counts = ones(L-W+1)
    for i in range(L-W):
        s = X[i:i+W]
        for j in range(1,W):
            t = X[i+j:i+j+W]
            if s == t:
                counts[i] = counts[i] + 1
                counts[i+j] = counts[i+j] + 1
    return counts


"""
The EM algorithm. 

Input:
Y, pygr database. dataset of sequences (As of 6/28/13, assume each sequence contains 1 subsequence
theta_motif, motif PWM matrix guess
theta_background_matrix, background PWM matrix guess
lambda_motif, motif frequency guess

Output:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency
"""
def Online_EM(Y, theta_motif, theta_background_matrix, lambda_motif):
    W = theta_motif.shape[0]#get the length of the motif
    #X = getSubsequences(Y,W)#this step may need to be removed for Online to save RAM
    #subsequences are grouped by sequences for normalization purposes
    #Is = [[sequenceToI(xij) for xij in xi] for xi in X]#list of indicator matrices, same dimensions as X
    s1_1 = lambda_motif#the expected number of occurrences of the motif
    s1_2 = theta_motif#the matrix holding the expected number of times a letter appears in each position, motif
    s2_2 = theta_background_matrix#the matrix holding the expected number of times a letter appears in each position, background
    n = 0#the counter
    nstart = 10000000#when to start averaging
    N = len(Y)#number of observations
    #reserve some memory. this was when each sequence had only one subsequence
    """fractions = zeros(N)
    pwms = zeros((N,W,4))
    backgrounds = zeros((N,4))
    """
    #prepare lists because we don't know how many subsequences we have in total
    fractions = list()
    pwms = list()
    backgrounds = list()
    expectations = list()
    fractions_sum = 0#should be deleted
    pwms_sum = zeros((W,4))#should be deleted
    backgrounds_sum = zeros((W,4))#should be deleted
    keys = Y.keys()
    shuffle(keys)
    #expectations.append(expected_LogLikelihood(Is, theta_motif, theta_background_matrix, lambda_motif))#add the expectation of the initial guess
    for y in keys:#iterate through each key in the FASTA file
        s = str(Y[y])#grab the whole sequence as a string
        L = len(s)#length of sequence
        starts = range(0,L-W+1)
        shuffle(starts)
        for start in starts:
            I = sequenceToI(s[start:start+W])#convert the subsequence to an indicator matrix
            step = 0.2*pow(n+1,-0.6)#the online step size. For OLO6a
            #step = 0.025*pow(n+1,-0.6)#the online step size. For OLO6a
            #step = 1.0/10000
            #E-step
            ds1_1 = Z0_I(I,theta_motif, theta_background_matrix,lambda_motif)
            #print y
            ds1_2 = ds1_1*I
            ds2_2 = (1-ds1_1)*I
            s1_1 = s1_1 + step*(ds1_1 - s1_1)
            s1_2 = s1_2 + step*(ds1_2 - s1_2)
            #print s1_2
            s2_2 = s2_2 + step*(ds2_2 - s2_2)
            #M-step
            lambda_motif = s1_1
            theta_motif = s1_2
            theta_motif = theta_motif/theta_motif.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob
            theta_background = s2_2.sum(axis = 0)#collapse the expected background counts into a single array
            theta_background = theta_background/theta_background.sum()#divide by the total counts to normalize to 1
            theta_background = array([theta_background])#prepare background for repeat
            theta_background_matrix = theta_background.repeat(W,axis=0)
            #save current parameters
            """fractions[n] = lambda_motif
            pwms[n] = theta_motif
            backgrounds[n] = theta_background
            fractions_sum = fractions_sum + lambda_motif
            pwms_sum = pwms_sum + theta_motif
            backgrounds_sum = backgrounds_sum + theta_background_matrix
            """
            fractions.append(lambda_motif)
            pwms.append(theta_motif)
            backgrounds.append(theta_background)
            #if n > nstart, then start using averaged parameters for the upcoming E-step
            #have to repeat the normalization to ensure probability is properly conserved
            if n > nstart:
                lambda_motif = mean(fractions[n/2:n],axis=0)#new fraction is mean of previous fractions
                theta_motif = mean(pwms[n/2:n],axis=0)#new pwm is mean of previous pwms
                theta_motif = theta_motif/theta_motif.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob
                theta_background = mean(backgrounds[n/2:n],axis=0)#new background is mean of previous backgrounds
                theta_background = theta_background/theta_background.sum()#divide by the total counts to normalize to 1
                #theta_background = array([theta_background])#prepare background for repeat
                theta_background_matrix = theta_background.repeat(W,axis=0)
                """
                lambda_motif = fractions_sum/(n+1)
                theta_motif = pwms_sum/(n+1)
                theta_motif = theta_motif/theta_motif.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob
                theta_background_matrix = backgrounds_sum/(n+1)
                theta_background_matrix = theta_background_matrix/theta_background_matrix.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob
                """
            #update the counter
            n = n + 1
            print n
            #the expected log likelihood, the objective function, based on current parameters
            #expectations.append(expected_LogLikelihood(Is, theta_motif, theta_background_matrix, lambda_motif))#add the expectation of the initial guess
    print s1_2
    #x = load('NRF1_Motif.npy')
    #pylab.plot([dist(x,y) for y in pwms])
    plot(fractions)
    show()
    #save('the_expectations', expectations)
    return theta_motif, theta_background_matrix, lambda_motif
    #E-step, this may be superfluous
    #Z, c0, c = E(I, theta_motif, theta_background_matrix, lambda_motif)
    #M-step, this may be superfluous
    #lambda_motif, theta_motif, theta_background_matrix = M(Z, n, c0, c)
    """
    for k in xrange(MAXITER):
        theta_motif_old = theta_motif#this is the only thing I need to save
        #E-step
        Z, c0, c = E(I, theta_motif, theta_background_matrix, lambda_motif)
        #M-step
        lambda_motif, theta_motif, theta_background_matrix = M(Z, n, c0, c)
        if dist(theta_motif, theta_motif_old) < TOL:
            break
    return lambda_motif, theta_motif, theta_background_matrix, k
    """

"""
The main online MEME algorithm.

Input:
Y, pygr database. dataset of sequences
W, width of motifs to search for
NPASSES, number of distinct motifs to search for

Assume that each sequence is the size of the motif for now.
That is, Y = X.
"""
def meme(Y,W,NPASSES):
    #6/28/13, check with initial conditions matching solution
    lambda_motif = 0.00625
    theta_motif = load('NRSF_test.npy')
    theta_uniform_background = array([[0.25, 0.25, 0.25, 0.25]])
    theta_uniform_background_matrix = theta_uniform_background.repeat(W,axis=0)#the initial guess for background is uniform distribution
    theta_motif, theta_background_matrix, lambda_motif = Online_EM(Y, theta_motif, theta_uniform_background_matrix, lambda_motif)
    outputMotif(lambda_motif, theta_motif, theta_background_matrix)
    
"""
Outputs the motif as a web logo.

Input:
lambda_motif - a double, fraction of subsequences 
theta_motif - a numpy array, the PWM
theta_background_matrix - a numpy array, the background model
"""
def outputMotif(lambda_motif, theta_motif, theta_background_matrix):
    print 'hello'
    print theta_motif
    data = LogoData.from_counts(counts=theta_motif,alphabet=unambiguous_dna_alphabet)
    options = LogoOptions()
    options.title = 'Motif'
    forma = LogoFormat(data, options)
    fout = open('results.png', 'w')
    png_formatter(data, forma, fout)
    fout.close()
    img = imread('results.png')
    imshow(img)
    show()
    #for now, just print, but will have to output a png later
    print lambda_motif
    print theta_motif
    print theta_background_matrix

if __name__ == "__main__":
    usage = "usage: %prog [options] <input FASTA>"
    description = "The program applies the Online MEME algorithm to find motifs in a FASTA file"
    parser = OptionParser(usage=usage,description=description)
    parser.add_option("-p", "--processes", help="optional number of parallelized processes")
    parser.add_option("-w", "--width", help="Width of the motif to search for", default="10")
    parser.add_option("-n", "--nummotifs", help="Number of motifs to search for", default="1")
    (options, args) = parser.parse_args()
    w = int(options.width)
    nmotifs = int(options.nummotifs)
    if len(args) == 1:#the program is correctly used, so do MEME
        sp = seqdb.SequenceFileDB(args[0])
        meme(sp,w,nmotifs)
        sp.close()
    else:
        parser.print_help()