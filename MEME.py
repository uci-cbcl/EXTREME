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
from numpy import load,inf, sign, dot, diag, array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, log, power, ceil, prod, zeros, concatenate
from numpy.random import rand
from itertools import repeat, chain, izip
from pygr import seqdb
from bisect import bisect_left


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
def pX_motif(X, theta_motif):
    d = {'A':0,'C':1,'G':2,'T':3}
    ps = [row[d[x]] for x, row in zip(X,theta_motif)]
    p = prod(ps)
    return p

"""
Equation 15 from the Bailey and Elkan paper. Calculates P(X|theta_background). That is,
the probability of the sequence given the background model. It is defined as the product
of the frequencies of the letter at each position.

Input:
X, A string. The sequence. Must only contain A, C, G, T
theta_background, Array of probabilities for each nucleotide, in order of A, C, G, T. Must sum to 1

Output:
p, A double. P(X|theta_background) 
"""
def pX_background(X, theta_background):
    d = {'A':0,'C':1,'G':2,'T':3}
    indices = [d[x] for x in X]
    ps = theta_background[indices]
    p = ps.prod()
    return p

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

Input:
X, a string. The sequence.
theta_motif, a matrix. The pwm of the motif.
theta_background, an array. The background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Notes:
5/6/13, work on this part
"""
def Z0(X,theta_motif, theta_background,lambda_motif):
    a = pX_motif(X,theta_motif)*lambda_motif#saves a calculation
    return a/(a + pX_background(X,theta_background)*(1-lambda_motif))


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
Equation 10 from the Bailey and Elkan paper. The expectation of the Z value for the sequence X,
given the motif frequency and the motif and background models. Note that g=2, where group 1
is the motif, and group 2 is the background. We can get the expectation of the Z value of 
the background by taking 1 - Z0.

Also computes a term in the summation, equation 11.

This is a modified version that uses the indicator matrix instead. Saves the trouble
of calculating the indicator matrix at each step. 

Input:
I, indicator matrix. Represents a sequence.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
Z0 - Expected value of Z for the the indicator matrix I
elogL - a summation term of equation 11
"""
def Z0_expLog_I(I,theta_motif, theta_background_matrix,lambda_motif):
    a = pI_motif(I,theta_motif)*lambda_motif#saves a calculation
    b = pI_background(I,theta_background_matrix)*(1-lambda_motif)#saves another calculation
    Z0 = a/(a + b)
    elogL = Z0*log(a) + (1-Z0)*log(b)
    return Z0,elogL

def expected_Log_Likelihood(Z,I):
    return 0

"""
Modified E-step of the EM algorithm. Accepts the list of indicator matrices I, the current motif
frequency lambda_motif, and both PWMs. Returns the expected Z values, and the 
expected number of times each letter appears c_jk.

Note that this uses the indicator matrices only. This modified version
also returns the expected log likelihood

Input:
I, a list of indicator matrices. This might get confusing with I being a list here and a matrix elsewhere.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
Z, a list of double. Same dimensions as I. Expected value of subsequence generated by motif
c0, an array. Equation 16. Expected number of times each letter appears generate by background
c, a matrix. Dimension Wx4. Equation 17. Expected number of times each letter appears at
each position generated by motif model.
expected_LogLikelihood, the expected log likelihood

"""
def Expectation(I, theta_motif, theta_background_matrix, lambda_motif):
    Z_and_elogL = [[Z0_expLog_I(Iij,theta_motif, theta_background_matrix,lambda_motif) for Iij in Ii] for Ii in I]
    Z = [[xij[0] for xij in xi] for xi in Z_and_elogL]
    #print Z
    expected_LogLikelihood = 0
    for xi in Z_and_elogL:
        for xij in xi:
            expected_LogLikelihood = expected_LogLikelihood + xij[1]
    L = theta_motif.shape[1]#Alphabet size same as number of columns in PWM
    c0 = zeros(L)
    c = zeros(theta_motif.shape)
    for Zi, Ii in izip(Z,I):
        for Zij, Iij in izip(Zi, Ii):
            c0 = c0 + (1 - Zij)*Iij.sum(axis=0)
            c = c + Zij*Iij#notice how this leaves the possibility of introducing erasing factors
    return Z, c0, c, expected_LogLikelihood

"""
The E-step of the EM algorithm. Accepts the list of indicator matrices I, the current motif
frequency lambda_motif, and both PWMs. Returns the expected Z values, and the 
expected number of times each letter appears c_jk.

Note that this uses the indicator matrices only. 

Input:
I, a list of indicator matrices. This might get confusing with I being a list here and a matrix elsewhere.
theta_motif, a matrix. The PWM of the motif.
theta_background_matrix, a matrix. Essentially a PWM of the background model.
lambda_motif, a double. The fraction of motifs among the sequences.

Output:
Z, a list of double. Same dimensions as I. Expected value of subsequence generated by motif
c0, an array. Equation 16. Expected number of times each letter appears generate by background
c, a matrix. Dimension Wx4. Equation 17. Expected number of times each letter appears at
each position generated by motif model.

"""
def E(I, theta_motif, theta_background_matrix, lambda_motif):
    Z = [[Z0_I(Iij,theta_motif, theta_background_matrix,lambda_motif) for Iij in Ii] for Ii in I]
    L = theta_motif.shape[1]#Alphabet size same as number of columns in PWM
    c0 = zeros(L)
    c = zeros(theta_motif.shape)
    for Zi, Ii in izip(Z,I):
        for Zij, Iij in izip(Zi, Ii):
            c0 = c0 + (1 - Zij)*Iij.sum(axis=0)
            c = c + Zij*Iij#notice how this leaves the possibility of introducing erasing factors
    return Z, c0, c

"""
The M-step of the EM algorithm. Accepts the expected Z values of each sequence and the expected
number of times each letter appears in each position (c_jk matrix), and returns the updated
motif frequency and PWM.

Input: 
Z, a list with the same dimensions as the subsequences list
n, the number of subsequences
c0, an array. Equation 16. Expected number of times each letter appears generate by background
c, a matrix. Dimension Wx4. Equation 17. Expected number of times each letter appears at
each position generated by motif model.

Output:
lambda_motif, a double. Motif frequency 
theta_motif, a matrix. PWM of the motif
theta_background_matrix, a matrix. Essentially a PWM of the background model
"""
def M(Z, n, c0, c):
    Z_total = 0
    for Zi in Z:
        for Zij in Zi:
            Z_total = Z_total + Zij
    lambda_motif = Z_total/n
    c0 = array([c0])
    c = concatenate((c0,c))
    f = dot(diag(1/c.sum(axis=1)),c)
    theta_motif = f[1:]
    theta_background = array([f[0]])
    theta_background_matrix = theta_background.repeat(theta_motif.shape[0],axis=0)
    return lambda_motif, theta_motif, theta_background_matrix

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
The EM algorithm. 

Input:
I, list of indicator matrices. Each matrix represents a subsequence
n, number of subsequences.
theta_motif, motif PWM matrix guess
theta_background_matrix, background PWM matrix guess
lambda_motif, motif frequency guess
TOL, tolerance for absolute change in PWM for convergence
MAXITER, maximum number of iterations

Output:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency
k, number of iterations performed
"""
def EM(I, n, theta_motif, theta_background_matrix, lambda_motif, TOL=1E-6, MAXITER=1000):
    #E-step, this may be superfluous
    #Z, c0, c = E(I, theta_motif, theta_background_matrix, lambda_motif)
    #M-step, this may be superfluous
    #lambda_motif, theta_motif, theta_background_matrix = M(Z, n, c0, c)
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
A troubleshoot version of the EM algorithm. This one prints the expected
log likelihood at each step. 

Input:
I, list of indicator matrices. Each matrix represents a subsequence
n, number of subsequences.
theta_motif, motif PWM matrix guess
theta_background_matrix, background PWM matrix guess
lambda_motif, motif frequency guess
TOL, tolerance for absolute change in PWM for convergence
MAXITER, maximum number of iterations

Output:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency
k, number of iterations performed
"""
def EM_troubleshoot(I, n, theta_motif, theta_background_matrix, lambda_motif, TOL=1E-6, MAXITER=1000):
    #E-step, this may be superfluous
    #Z, c0, c = E(I, theta_motif, theta_background_matrix, lambda_motif)
    #M-step, this may be superfluous
    #lambda_motif, theta_motif, theta_background_matrix = M(Z, n, c0, c)
    for k in xrange(MAXITER):
        theta_motif_old = theta_motif#this is the only thing I need to save
        #E-step
        Z, c0, c, expected_LogLikelihood = Expectation(I, theta_motif, theta_background_matrix, lambda_motif)
        print "Step " + str(k) + ":" + " " + str(expected_LogLikelihood)
        print theta_motif
        #M-step
        lambda_motif, theta_motif, theta_background_matrix = M(Z, n, c0, c)
        if dist(theta_motif, theta_motif_old) < TOL:
            break
    return lambda_motif, theta_motif, theta_background_matrix, k

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
Samples Q (ceiled) random subsequences from I, the list of indicator matrices.
This function is needed because I is a ragged array and therefore cannot
be directly manipulated.

Input:
I, list of indicator matrices to be sampled from
nlist, cumulative summed index ends of each sequence
n, number of subsequences
Q, number of subsequences to sample. Less than n

Output:
Is, list of Q sampled indicator matrices
"""
def sampleIndicator(I, nlist, n, Q):
    indices = sample(xrange(n),Q)#sampled indices
    seqindices = [bisect_left(nlist, ind) for ind in indices]
    Is = [I[seqind][ind - nlist[seqind] - 1] for seqind, ind in zip(seqindices,indices)]#which sequences to use
    #note the subtraction. we are indexing from the right here
    return Is

"""
The relative entropy per column, RE/col. A measure of the "crispness" of the motif.

Input:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix

Output:
recol, the relative entropy per column
"""
def relEntCol(theta_motif, theta_background_matrix):
    W = theta_motif.shape[0]
    recol = 1.0/W*(theta_motif*log(theta_motif/theta_background_matrix)).sum()
    return recol

"""
Equation 47. Used for evaluation in the bisection method to find the roots of
the function.

Input:
m - the value to evaluate the function at
gamma - the user-defined fraction of entropy. default:0.1

Output:
Evaluated function
"""
def objFunction(m,gamma):
    L = 4.0
    p = 1/L
    return (1-m)*log((1-m)/(L-1)) - log(p) + m*log(m) + gamma*log(p)

"""
Bisection method to solve the equation for m.

Input:
gamma - used defined fraction of entropy. default:0.1
a - initial guess for left boundary
c - initial guess for right boundary
tol - tolerance for difference between boundaries
nmax - max number of iterations

Output:
b - the root
"""
def bisection(gamma,a,c,tol=10^-6,nmax=100):
    n = 0
    while n < nmax:
        b = (a+c)/2
        f = objFunction(b,gamma)
        if (f == 0 or (c-a)/2 < tol):
            break
        if sign(f) == sign(objFunction(a,gamma)):
            a = b
        else:
            c = b
        n = n + 1
    return b



"""
Equations 45-47. Finds the initial position weight matrix to use as a starting point, given a
sampled subsequence. Because the value of m is constant reused, it is pre-calculated.

Input:
I - the sampled indicator matrix to be converted to a starting point
m - the value derived from the bisection method

Output:
theta_motif_guess - the starting point derived from I
"""
def startingPoint(I, m):
    a = m*I
    b = (1.0-m)/(4.0-1)*(I==False)
    theta_motif_guess = a+b
    return theta_motif_guess
    
"""
Outputs the motif as a web logo.

Input:
lambda_motif - a double, fraction of subsequences 
theta_motif - a numpy array, the PWM
theta_background_matrix - a numpy array, the background model
"""
def outputMotif(lambda_motif, theta_motif, theta_background_matrix):
    c = theta_motif.T
    d = {'A':c[0],'C':c[1],'G':c[2],'T':c[3]}
    m = motifs.Motif(alphabet=IUPAC.unambiguous_dna,counts=d)
    b = theta_background_matrix[0]
    back = {'A':b[0],'C':b[1],'G':b[2],'T':b[3]}
    m.background = back
    m.weblogo('results.png')
    #for now, just print, but will have to output a png later
    print lambda_motif
    print theta_motif
    print theta_background_matrix

"""
The main MEME algorithm.

Input:
Y, pygr database. dataset of sequences
W, width of motifs to search for
NPASSES, number of distinct motifs to search for

5/14/13:
Assume that each sequence is the size of the motif for now.
That is, Y = X.
"""
def meme(Y,W,NPASSES):
    X = getSubsequences(Y,W)#this step may need to be moved for Online to save RAM
    #subsequences are grouped by sequences for normalization purposes
    I = [[sequenceToI(xij) for xij in xi] for xi in X]#list of indicator matrices, same dimensions as X
    nlist = array([len(x) for x in X])#a list of the number of subsequences corresponding to each sequence
    n = nlist.sum()#number of subsequences. this step may need to be changed for Online
    print n
    nlist = nlist.cumsum()#now nlist can be searched by bisection
    nlist = nlist - 1#the subtraction is because I need indices, not lengths, for the random sampling
    N = len(Y)#number of sequences
    alpha = 0.99#probability of catching sequence with motif
    m = bisection(0.1,0.25,1)#m, for initial PWM guesses. The boundaries are determined from the minimum of the function
    theta_uniform_background = array([[0.25, 0.25, 0.25, 0.25]])
    theta_uniform_background_matrix = theta_uniform_background.repeat(W,axis=0)#the initial guess for background is uniform distribution
    """
    print 'Running EM'
    lambda_motif = 0.01
    theta_motif = load('NRF1_Motif.npy')
    theta_motif[0] = array([0.25,0.25,0.25,0.25])
    lambda_motif, theta_motif, theta_background_matrix, k = EM_troubleshoot(I, n, theta_motif, theta_uniform_background_matrix, lambda_motif)
    print 'Iterations'
    print k
    outputMotif(lambda_motif, theta_motif, theta_background_matrix)
    """
    for i in range(NPASSES):#search for NPASSES of motifs of width W
        a = arange(0,log2(n/(2*W*sqrt(N)))+1,1)#the +1 ensures the last guess included
        b = a + log2(sqrt(N)/n) 
        lambda0s = power(2,b)#array holding the heuristic initial lambda values to try
        lambda0s = lambda0s[-1:]#for now only want the last guessed lambda
        lambda0s = [0.3]
        print lambda0s#troubleshooting
        Qs = [int(ceil(log2(1-alpha)/log2(1-lambda0))) for lambda0 in lambda0s]
        print Qs#for troubleshooting
        best_expectedLog = -inf
        #this loop and inner loop finds the best starting point
        for lambda0 in lambda0s:
            Q = int(ceil(log2(1-alpha)/log2(1-lambda0)))#number of sequences to sample to catch motif with prob alpha
            print 'Q'
            print Q
            print 'lambda'
            print lambda0
            sampledIndicators = sampleIndicator(I, nlist, n, Q)#guarantee to sample at least one sequence with motif with alpha confidence
            dx = 1
            for sampledI in sampledIndicators:
                #derive a PWM guess from the sampled subsequence
                theta_motif_guess = startingPoint(sampledI,m)
                #E-step
                Z, c0, c = E(I, theta_motif_guess, theta_uniform_background_matrix, lambda0)
                #M-step
                lambda_motif_updated, theta_motif_updated, theta_background_matrix_updated = M(Z, n, c0, c)
                #Completed 1 iteration, now do E-step to get the likelihood from 1 iteration
                Z, c0, c, expected_logLikelihood = Expectation(I, theta_motif_updated, theta_background_matrix_updated, lambda_motif_updated)
                #pick the best guess for the parameters, which
                #will be used for the full EM iterations
                if expected_logLikelihood > best_expectedLog:
                    best_expectedLog = expected_logLikelihood
                    print expected_logLikelihood
                    lambda_motif = lambda_motif_updated#current best guess fraction
                    theta_motif = theta_motif_updated#current best guess PWM
                    theta_background_matrix = theta_background_matrix_updated#current best guess background model
                print dx
                dx = dx + 1
            #now run EM algorithm on the current best guess
            print 'Running EM'
            lambda_motif, theta_motif, theta_background_matrix, k = EM_troubleshoot(I, n, theta_motif, theta_background_matrix, lambda_motif)
            print 'Iterations'
            print k
            outputMotif(lambda_motif, theta_motif, theta_background_matrix)
            
                
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
    if len(args) == 1:#the program is correctly used, so do MEME
        sp = seqdb.SequenceFileDB(args[0])
        meme(sp,10,1)
        sp.close()
    else:
        parser.print_help()