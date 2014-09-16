#!/usr/bin/env python
from argparse import ArgumentParser
import random
import copy
import errno
import sys
import sequence
from collections import deque
from numpy import round_,mean,load,save,inf, sign, dot, diag, array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, log, power, floor, ceil, prod, zeros, ones, concatenate, argmin
from itertools import chain

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
Z0 - Expected value of Z for the the indicator matrix I. Returns 0 if I is None
"""
def Z0_I(I,theta_motif, theta_background_matrix,lambda_motif, fudgefactor=1.0):
    if I is None:
        return 0
    a = fudgefactor*pI_motif(I,theta_motif)*lambda_motif#saves a calculation
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
Each column represents a nucleotide, in alphabetical order. If X contains a deleted base pair,
the returned indicator matrix will be a None.
"""
def sequenceToI(X):
    if 'N' in X:
        return None
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
    return [[y[k:k+W] for k in xrange(len(y)-W+1)] for y in Y]



"""
Calculates and returns the symmetrized KL Divergence between two PFMs.
"""
def KLD(x,y):
    Z = 0.5*(sum(x*log(x/y)) + sum(y*log(y/x)))
    return Z

def get_probs(seqs, alphabet_string):
    """ Get the observed probabilities of the letters in a set
    of sequences.  Ambiguous characters are ignored.
    Uses an "add-one" prior. From dreme.py by T. Bailey"""

    freqs = {}
    # initialize with add-one count
    for char in alphabet_string:
        freqs[char] = 1
    # get the frequencies of DNA letters in the sequences
    for seq in seqs:
        for char in seq:
            if freqs.has_key(char):
                freqs[char] += 1
            else:
                freqs[char] = 0         # ambiguous letter
    # get the total number of non-ambiguous letters
    n = 0.0
    for char in alphabet_string:
        n += freqs[char]
    # normalize the probabilities
    probs = {}
    for char in alphabet_string:
        probs[char] = freqs[char]/n

    return probs



"""
The online EM algorithm. 

Input:
Is, list of lists of indicator matrices. dataset of sequences
seqindpairs, list of tuples of valid sequence and subsequence indices to use
theta_motif, motif PWM matrix guess
theta_background_matrix, background PWM matrix guess
lambda_motif, motif frequency guess
smoothing, whether to smooth (default: False)
revcomp, whether to use both strands (default: True)
B, pseudo-counts parameter (default: 0.001)

Output:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency
"""
def Online_EM(Is, seqindpairs, theta_motif, theta_background_matrix, lambda_motif, fudgefactor, minsites, maxsites, initialstep=0.05, B=0.0001, smoothing=False, revcomp=True):
    W = theta_motif.shape[0]#get the length of the motif
    s1_1 = lambda_motif#the expected number of occurrences of the motif
    s1_2 = theta_motif#the matrix holding the expected number of times a letter appears in each position, motif
    s2_2 = theta_background_matrix#the matrix holding the expected number of times a letter appears in each position, background
    n = 0#the counter
    N = len(Is)#number of sequences
    pwms = list()
    backgrounds = list()
    expectations = list()
    fractions_sum = 0#should be deleted
    pwms_sum = zeros((W,4))#should be deleted
    backgrounds_sum = zeros((W,4))#should be deleted
    firstmotif = theta_motif#for tracking how much the motif changes
    mu = theta_background_matrix#the first background matrix is the average frequencies in the negative set
    Bmu = B*mu#the priors to be added each step
    g0 = max(initialstep,lambda_motif*10)
    g1 = -0.6
    print "Initial step size of " + str(g0)
    print "Running Online EM algorithm..."
    pwm_deque = deque(maxlen=100)
    for ps in range(5):
        for seqindpair in seqindpairs:#iterate through each sequence index and start pair
            seqind = seqindpair[0]
            start = seqindpair[1]
            step = g0*pow(n+1,g1)#the online step size. For OLO6a
            I = Is[seqind][start]#grab the current indicator matrix
            Z = Z0_I(I,theta_motif, theta_background_matrix,lambda_motif, fudgefactor)#perhaps implement RC better here?
            this_strand = True
            if revcomp:#if the user wants reverse complements
                Ir = I_rc(I)
                Zr = Z0_I(Ir,theta_motif, theta_background_matrix,lambda_motif, fudgefactor)
                if Zr > Z:
                    Z = Zr
                    I = Ir
                    this_strand = False#opposite strand strong, use reverse complement
            ds1_1 = Z
            I = I + Bmu
            ds1_2 = ds1_1*I
            ds2_2 = (1-ds1_1)*I
            s1_1 = s1_1 + step*(ds1_1 - s1_1)
            s1_2 = s1_2 + step*(ds1_2 - s1_2)
                #print s1_2
            s2_2 = s2_2 + step*(ds2_2 - s2_2)
                #M-step
            lambda_motif = s1_1
            theta_motif = s1_2# + step*Bmu
            theta_motif = theta_motif/theta_motif.sum(axis=1)[:,newaxis]#ensures each row has sum 1, for prob

            theta_background = s2_2.sum(axis = 0)#collapse the expected background counts into a single array
            theta_background = theta_background/theta_background.sum()#divide by the total counts to normalize to 1
            theta_background = array([theta_background])#prepare background for repeat
            theta_background_matrix = theta_background.repeat(W,axis=0)

            #update the counter
            n = n + 1
            pwm_deque.append(theta_motif)
            #print n
            #the expected log likelihood, the objective function, based on current parameters
            #expectations.append(expected_LogLikelihood(Is, theta_motif, theta_background_matrix, lambda_motif))#add the expectation of the initial guess
        print "KLD:",KLD(pwm_deque[0],pwm_deque[-1])
        if KLD(pwm_deque[0],pwm_deque[-1]) < 1e-6:
            print "KLD threshold met on pass",ps+1
            break
        else:
            print "KLD threshold not met. Doing another pass"
            g1 = (g1-1)/2
            #if lambda_motif*N < minsites or lambda_motif*N > maxsites:
            #    print lambda_motif*N
            #    print "But it probably was a bad run anyways, so break loop"
                #break
    return theta_motif, theta_background_matrix, lambda_motif



"""
A modified version of the extreme algorithm. It uses PFMs from a database such as JASPAR
or the UW ENCODE group, and the predicted number of sites, as seeds.

Input:
Y, list of strings. dataset of sequences
minsites, the minimium number of sites
maxsites, the maximum number of sites. If 0, it is automatically changed to 5 times the number of predicted sites
pwm_guess, the PFM of the initial guess
tries, number of different "fudge factors"/bias factors to try before giving up
Output:
fractions
"""
def extreme(Y,neg_seqs,minsites,maxsites,pwm_guess,initialstep=0.05,tries=15,revcomp=True):
    #6/28/13, check with initial conditions matching solution
    #p = Pool(64)
    #s=p.map(functools.partial(f,y=Y),range(64))
    #print s
    BIGLOG = log(sys.float_info.max)
    pos_seqs = Y
    _dna_alphabet = 'ACGT'
    #The discovered motifs and results to be reported
    discovered_lambda_motifs = list()
    discovered_theta_motifs = list()
    discovered_theta_background_matrices = list()
    discovered_logevs = list()
    discovered_nonoverlapsites = list()
    #All motifs and results to be reported
    all_lambda_motifs = list()
    all_theta_motifs = list()
    all_theta_background_matrices = list()
    all_logevs = list()
    dprobs = get_probs(neg_seqs, _dna_alphabet)
    #generate first order Markov background based on nucleotide frequencies
    print 'Getting background model'
    theta_background = array([[dprobs['A'], dprobs['C'], dprobs['G'], dprobs['T']]])
    #print theta_background
    #lists to hold the motifs and results in this round
    lambda_motifs = list()
    theta_motifs = list()
    theta_background_matrices = list()
    fractionses = list()
    distanceses = list()
    logevs = list()
    DQ = pwm_guess
    W = DQ.shape[0]
    #print 'Using starting point from DREME PWM generation...'
    #n = sum([max(0,len(y) - W + 1) for y in Y])#gets number of subsequences
    print 'Getting subsequences'
    X = getSubsequences(Y,W)#this step may need to be removed for Online to save RAM
    #subsequences are grouped by sequences for normalization purposes
    print 'Getting indicator matrices'
    Is = [[sequenceToI(xij) for xij in xi] for xi in X]#list of indicator matrices for this specific W, same dimensions as X    
    DR = theta_background.repeat(DQ.shape[0],axis=0)#the initial guess for background is uniform distribution
    print "Scanning sequence with current PWM guess"
    pos = guess_positive_sites(DQ, DR, Is)
    print "Guessing",pos,"sites"
    #print "Found",pos,"consensus sequence matches in the positive sequences"
    #print "Found",neg,"consensus sequence matches in the negative sequences"
    if maxsites == 0:
        maxsites = pos*5
        print "Maximum number of sites not specified, so setting it to",maxsites
    #Valid sequence and subsequence index combinations to search. 
    #subsequences with deleted base pairs will yield a None indicator matrix, which will be ignored in this list
    seqindpairs = [[(seqind,subseqind) for subseqind in xrange(len(Is[seqind])) if Is[seqind][subseqind] is not None] for seqind in xrange(len(Is))]
    seqindpairs = list(chain(*seqindpairs))
    n = len(seqindpairs)#total number of subsequences
    random.shuffle(seqindpairs)
    #The bounds for the fudge factor
    a = 0.0
    b = 1.0
    c = 1.0
    for t in range(tries):
        print 'Try ' + str(t + 1)
        print 'Using a fudge factor of ' + str(b)
        fudgefactor = b
        theta_motif = DQ
        lambda_motif = 1.0*pos/n#guess twice the number regular expression matches
        theta_background_matrix = theta_background.repeat(theta_motif.shape[0],axis=0)#the initial guess for background is uniform distribution
        theta_motif, theta_background_matrix, lambda_motif = Online_EM(Is, seqindpairs, theta_motif, theta_background_matrix, lambda_motif, fudgefactor, minsites, maxsites, initialstep)
        print 'Finding number of motif sites'
        if lambda_motif < 1e-9:
            nsites_dis = 0
        else:
            nsites_dis = get_nsites_dis(theta_motif, theta_background_matrix, lambda_motif, Is)
        print 'Found ' + str(nsites_dis) + ' sites'
        #if there are too many discovered sites, something is wrong, so assign a high E-value
        shouldIBreak = False
        if nsites_dis > maxsites:#for now, assume problem if more than 10 instances per sequence
            print 'Too many sites found. Setting log E-value to max value. Lowering fudge factor and reshuffling'
            logev = BIGLOG
            c = b
            b = mean([a,b])
            random.shuffle(seqindpairs)
        elif nsites_dis < minsites:
            print 'Not enough sites found. Setting log E-value to max value. Raising fudge factor and reshuffling'
            logev = BIGLOG
            a = b
            b = mean([b,c])
            random.shuffle(seqindpairs)
        else:
            shouldIBreak = True
            print 'Motif has an acceptable number of sites'    
            try:
                import meme as me
                mm = me.MEME(theta_motif, theta_background_matrix[0], lambda_motif, Y, nsites_dis)
                print 'Calculating log E-value'
                mm.calc_ent()
                logev = mm.get_logev()
                print 'Log E-value: ' + str(logev)
            except ImportError:
                print "You did not install the MEME Cython wrapper, so the E-value will be set to the largest possible value"
                logev = BIGLOG
            logevs.append(logev)
            lambda_motifs.append(lambda_motif)
            theta_motifs.append(theta_motif)
            theta_background_matrices.append(theta_background_matrix)
        #save everything
        all_logevs.append(logev)
        all_lambda_motifs.append(lambda_motif)
        all_theta_motifs.append(theta_motif)
        all_theta_background_matrices.append(theta_background_matrix)
        if shouldIBreak:
            break
    #went through all tries or found a motif with acceptable number of sites
    #if no valid motif found, then exit
    if len(logevs) == 0:
        print 'No motif found, so do nothing...'
    best_index = argmin(logevs)
    best_theta_motif = theta_motifs[best_index]
    best_theta_background_matrix = theta_background_matrices[best_index]
    best_lambda_motif = lambda_motifs[best_index]
    best_logev = logevs[best_index]
    print "Using PWM with log EV: " + str(best_logev)
    discovered_lambda_motifs.append(best_lambda_motif)
    discovered_theta_motifs.append(best_theta_motif)
    discovered_theta_background_matrices.append(best_theta_background_matrix)
    discovered_logevs.append(best_logev)
    pos_nsites, neg_nsites = erase_motif(best_theta_motif, best_theta_background_matrix, best_lambda_motif, pos_seqs, neg_seqs)
    discovered_nonoverlapsites.append(pos_nsites)
    return all_theta_motifs, all_theta_background_matrices, all_lambda_motifs, all_logevs, \
        discovered_theta_motifs, discovered_logevs, discovered_nonoverlapsites


"""
For guesssing the number of motif matches for the initial lambda_m. Calculates
goodness-of-fit score, G, for each subsequence (or rather, indicator matrix)
and checks whether G exceeds the threshold (default: 0.7).

Input: 
theta_motif, the PWM (assumed to be trimmed already)
theta_background_matrix, background frequencies, same size as theta_motif
Is, list of list of indicator matrices

Output:
pos, guess for the number of positive sites
"""
def guess_positive_sites(theta_motif, theta_background_matrix, Is, Gthresh=0.7):
    logodds_matrix = log(theta_motif/theta_background_matrix)#spec matrix
    Vmax = logodds_matrix.max(axis=1).sum()
    pos = 0
    for J in Is:
        for I in J:
            if I is None:
                continue
            G = max(goodness_fit(logodds_matrix, Vmax, I),goodness_fit(logodds_matrix, Vmax, I_rc(I)))
            pos += G > Gthresh
    return pos


"""
Calculates goodness of fit, G.

Input:
logodds_matrix, the log odds matrix
I, an indicator matrix for a subsequence
Vmax, maximum log-odd score

Output:
G - the goodness of fit
"""
def goodness_fit(logodds_matrix, Vmax, I):
    V = sum(logodds_matrix*I)
    G = V/Vmax
    return G

"""
Searches for motif instances in the positive and negative sequences and replaces
with N's. It should be noted that as the motif has been trimmed, the score will
be reduced (scores are purely additive by column), and the threshold is not
adjusted for the trimming.

Input: 
theta_motif, the PWM (assumed to be trimmed already)
theta_background_matrix, background frequencies, same size as theta_motif
lambda_motif, fraction of subsequences that are generated by motif
pos_seqs, list of positive sequences
neg_seqs, list of negative sequences

Output:
Updated positive and negative sequences with motif sites deleted
Also output number of sites erased in both sequence sets
"""
def erase_motif(theta_motif, theta_background_matrix, lambda_motif, pos_seqs, neg_seqs, revcomp=True):
    t = log((1-lambda_motif)/lambda_motif)#Threshold
    spec = log(theta_motif/theta_background_matrix)#spec matrix
    W = theta_motif.shape[0]#width of the motif
    ens = W * 'N'#N sequence to replace motif sites with
    pos_nsites_dis = 0
    print 'Erasing motif from positive sequences'
    for i in range(len(pos_seqs)):
        s = pos_seqs[i]#grab the string
        L = len(s)
        end = L - W + 1#go to last W-mer of sequence
        j = 0
        while j < end:
            subs = s[j:j+W]#grab the subsequence
            if 'N' in subs:#ignore subsequences with deleted portions
                j += 1
                continue
            I = sequenceToI(subs)#subsequence has no deleted portions, so make indicator matrix
            if revcomp:
                #a is boolean that tells whether a hit is found
                a = sum(spec*I) > t or sum(spec*I_rc(I)) > t
            else:
                a = sum(spec*I) > t
            if a:#hit found, increment, erase, and move index
                pos_nsites_dis += 1
                s = s[0:j] + ens + s[j+W:]
                j += W
            else:#no hit found, move index up by 1
                j += 1
        pos_seqs[i] = s#done erasing, store result
    print 'Erased ' + str(pos_nsites_dis) + ' sites from the positive sequences'     
    print 'Erasing motif from negative sequences'   
    neg_nsites_dis = 0
    for i in range(len(neg_seqs)):
        s = neg_seqs[i]#grab the string
        L = len(s)
        end = L - W + 1#go to last W-mer of sequence
        j = 0
        while j < end:
            subs = s[j:j+W]#grab the subsequence
            if 'N' in subs:#ignore subsequences with deleted portions
                j += 1
                continue
            I = sequenceToI(subs)#subsequence has no deleted portions, so make indicator matrix
            if revcomp:
                #a is boolean that tells whether a hit is found
                a = sum(spec*I) > t or sum(spec*I_rc(I)) > t
            else:
                a = sum(spec*I) > t
            if a:#hit found, increment, erase, and move index
                neg_nsites_dis += 1
                s = s[0:j] + ens + s[j+W:]
                j += W
            else:#no hit found, move index up by 1
                j += 1
        neg_seqs[i] = s#done erasing, store result
    print 'Erased ' + str(neg_nsites_dis) + ' sites from the negative sequences'
    return (pos_nsites_dis, neg_nsites_dis)


"""
Given an indicator matrix, return its reverse complement.

Input:
I, indicator matrix

Output:
I_rc, reverse complement of the indicator matrix
"""
def I_rc(I):
    return I[::-1,::-1]	

"""
Finds the number of discrete motif sites. Uses the threshold decision theory strategy
from Bailey and Elkan. Currently searches both strands.

Input:
Is, list of lists of indicator matrices. dataset of sequences
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency

Output:
nsites_dis, integer number of discovered motif sites

"""
def get_nsites_dis(theta_motif, theta_background_matrix, lambda_motif, Is, revcomp=True):
    t = log((1-lambda_motif)/lambda_motif)#Threshold
    spec = log(theta_motif/theta_background_matrix)#spec matrix
    nsites_dis = 0#discrete sites discovered
    for s in Is:
        for I in s:
            if I is None:#if this location intersects with a deleted base pair, go to next loop step
                continue
            if revcomp:
                a = sum(spec*I) > t or sum(spec*I_rc(I)) > t
            else:
                a = sum(spec*I) > t
            nsites_dis += a
    return nsites_dis

# print very large or small numbers
# from dreme.py by T. Bailey
def sprint_logx(logx, prec, format):
    """ Print x with given format given logx.  Handles very large
    and small numbers with prec digits after the decimal.
    Returns the string to print."""
    log10x = logx/log(10)
    e = floor(log10x)
    m = pow(10, (log10x - e))
    if ( m + (.5*pow(10,-prec)) >= 10):
        m = 1
        e += 1
    str = format % (m, e)
    return str

"""
Outputs the motif as a web logo. Saves fractions as .npy.

Input:
lambda_motif - a double, fraction of subsequences 
theta_motif - a numpy array, the PWM
theta_background_matrix - a numpy array, the background model
logev - log E-value of motif
k - the motif number index
outstr - the prefix for the output files
"""
def outputMotif(theta_motif, theta_background_matrix, lambda_motif, logev, k, outstr):
    from weblogolib import LogoData, LogoOptions, LogoFormat, png_formatter, eps_formatter, unambiguous_dna_alphabet
    _pv_format = "%3.1fe%+04.0f"
    f_string = sprint_logx(log(lambda_motif), 1, _pv_format)
    g_string = sprint_logx(logev, 1, _pv_format)
    print(("Motif {0:s} had a fraction of {1:s}").format(str(k), f_string))
    print(("Motif {0:s} had an E-value of {1:s}").format(str(k), g_string))
    print 'Saving motif as a png...'
    data = LogoData.from_counts(counts=theta_motif,alphabet=unambiguous_dna_alphabet)#,prior=theta_background_matrix[0])#Does prior mess things up?
    options = LogoOptions()
    options.title = 'Motif'
    forma = LogoFormat(data, options)
    fout = open(outstr + "Motif_" + str(k) + '.png', 'w')
    png_formatter(data, forma, fout)
    fout.close()
    print 'Saving motif as an eps...'
    fout = open(outstr + "Motif_" + str(k) + '.eps', 'w')
    eps_formatter(data, forma, fout)
    fout.close()

    
"""
Outputs the discovered motifs in the MEME format

Input:
disc_pwms - list of numpy arrays of the PWMs
disc_logevs - list of the log E-valus of the discovered motifs
disc_nsites - list of motif sites
outpre - the prefix of the output files
"""
def outputMEMEformat(disc_pwms, disc_logevs, disc_nsites, outpre):
    _pv_format = "%3.1fe%+04.0f"
    f = open(outpre+"MEMEoutput.meme","w")
    f.write("MEME version 4.9.0\n\n")
    f.write("ALPHABET= ACGT\n\n")
    f.write("strands: + -\n\n")
    f.write("Background letter frequencies (from uniform background):\n")
    f.write("A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n")
    n = 1
    for pwm, logev, nsites in zip(disc_pwms, disc_logevs, disc_nsites):
        g_string = sprint_logx(logev, 1, _pv_format)
        x = round_(pwm,3)
        w = x.shape[0]
        y = str(x).replace('[','').replace(']','').replace('  ',' ').replace('1.  ','1.000').replace('0.  ','0.000').replace('0.\n', '0.000\n').replace('1.\n', '1.000\n').replace('\n  ','\n')[1:]
        if y[-1] == '.':
            y += '000'#add zeros at end if there's a hanging non-decimal number
        f.write('MOTIF M' + str(n) + ' O'+str(n)+'\n\n')
        f.write('letter-probability matrix: alength= 4 w= ' + str(w) + ' nsites= ' + str(nsites) + ' E= ' + g_string +  '\n')
        f.write(' ' + y)
        f.write('\n\n')
        n += 1
    f.close()        

"""
The main executable function
"""
def main():
    usage = "usage: %prog [options] <input FASTA>"
    description = "The program applies a modified EXTREME algorithm to find motifs in a FASTA file. It accepts a positive sequence set, a negative sequence set, a list of seed PFMs, and an index number indicating which of the seed PFMs to use"
    parser = ArgumentParser(description=description)
    parser.add_argument('fastafile', metavar='f', help='FASTA file containing the sequences')
    parser.add_argument('negfastafile', metavar='g', help='Negative FASTA file. This is for comparison so that you know the motif you discovered is over-represented.')
    parser.add_argument('jfile', metavar='j', help='File containing PWM seeds')
    parser.add_argument('indexvalue', metavar='i', help='Which seed from the Minimal MEME Format file to use (it is an integer ranging from 1 to the total number of PFM seeds in your file)', type=int)
    parser.add_argument("-p", "--pseudocounts", help="Pseudo counts added to initial PFM guess. Default:0.0", type=float, default=0.0)
    parser.add_argument("-q", "--initialstep", help="The initial step size for the online EM algorithm. A VERY sensitive parameter. I get best success for ChIP size data (about 100,000 to 1,000,000 bps) with a step size of 0.05. For DNase footprinting, which usually has >5,000,000 bps, I find 0.02 works best. Default:0.05", type=float, default=0.05)    
    parser.add_argument("-maxsites", dest="maxsites", help="Maximum number of expected sites for the motif. If not specified, defaults to 5 times number of initial predicted sites.", type=int, default=0)
    parser.add_argument("-minsites", dest="minsites", help="Minimum number of expected sites for the motif. Default: 10", type=int, default=10)
    parser.add_argument("-t", "--tries", dest="tries", help="Number of tries for each motif discovered. The fudge factor is changed until the number of discovered sites is in the \"acceptable\" range", type=int, default=15)
    parser.add_argument("-s", "--seed", dest="seed", help="Random seed", type=int, default=1)
    parser.add_argument("-saveseqs", "--saveseqs", dest="saveseqs", help="If specified, save sequences to current directory", action='store_true')
    import time
    print "Started at:"
    print time.ctime()
    starttime = time.time()
    args = parser.parse_args()
    seed = args.seed
    initialstep = args.initialstep
    minsites = args.minsites
    maxsites = args.maxsites
    random.seed(seed)
    jfile = open(args.jfile,'r')
    from numpy import fromstring
    from string import join
    lines = jfile.readlines()
    j = 0
    for i in range(len(lines)):
        line = lines[i]
        if '>' in line:#This is a name line, so read in next lines for matrix
            j += 1
            if j == args.indexvalue:#at the desired index
                parts = lines[i].split()
                pos_cs = parts[1]
                motifname = parts[0][1:]
                w = len(pos_cs)
                strlines = lines[i+1:i+1+w]
                pwm_string = ''
                for strline in strlines:
                    strparts = strline.split()
                    for strpart in strparts:
                        pwm_string += strpart + ' '
                #print pwm_string
                pwm_guess = fromstring(pwm_string,sep=' ',dtype=float)
                pwm_guess = pwm_guess.reshape((w,4))
                break
    print 'Using initial motif guess',motifname
    print 'Adding',str(args.pseudocounts),'pseudocounts and normalizing'
    pwm_guess = pwm_guess + args.pseudocounts
    pwm_guess = pwm_guess/pwm_guess.sum(axis=1)[:,newaxis]
    jfile.close() 
   
    # make the directory (recursively)
    import os
    outdir = motifname
    outpre = outdir + "/"
    clobber = True
    try:#adapted from DREME.py by T. Bailey
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not clobber:
                print >> sys.stderr, ("output directory (%s) already exists "
                "but EXTREME was not told to clobber it") % (outdir); sys.exit(1)
        else: raise
    #Use DREME's SeqIO to read in FASTA to list
    seqs = sequence.convert_ambigs(sequence.readFASTA(args.fastafile, None, True))
    #print seqs
    negseqs = sequence.convert_ambigs(sequence.readFASTA(args.negfastafile, None, True))
    tries = args.tries
    theta_motifs, theta_background_matrices, lambda_motifs, logevs, disc_pwms, disc_logevs, disc_nsites = extreme(seqs,negseqs,minsites,maxsites,pwm_guess,initialstep,tries)
    k = 1
    outputMEMEformat(disc_pwms, disc_logevs, disc_nsites, outpre)
    try:
        from weblogolib import LogoData, LogoOptions, LogoFormat, png_formatter, eps_formatter, unambiguous_dna_alphabet
        for theta_motif, theta_background_matrix, lambda_motif, logev in zip(theta_motifs, theta_background_matrices, lambda_motifs, logevs):
            outputMotif(theta_motif, theta_background_matrix, lambda_motif, logev, k, outpre)
            k = k+1
    except ImportError:
        print "You do not have Weblogolib, so sequence logos will not be made"
    
    
    if args.saveseqs:
        print "Saving Positive sequences to Positive_seq.fa"
        pos_file = open("Positive_seq.fa","w")
        for s in range(len(seqs)):
            pos_file.write(">sequence"+str(s+1)+"\n")
            pos_file.write(seqs[s]+"\n")
        pos_file.close()
        print "Saving Negative sequences to Negative_seq.fa"
        neg_file = open("Negative_seq.fa","w")
        for s in range(len(negseqs)):
            neg_file.write(">sequence"+str(s+1)+"\n")
            neg_file.write(negseqs[s]+"\n")
        neg_file.close()
    print "Ended at:"
    print time.ctime()
    stoptime = time.time()
    duration = stoptime - starttime
    print "Duration:", duration


if __name__=='__main__':
    main()
