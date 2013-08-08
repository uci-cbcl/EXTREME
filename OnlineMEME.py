#!/usr/bin/env python
import os
from argparse import ArgumentParser
from multiprocessing import Pool
from numpy import mean,load,save,inf, sign, dot, diag, array, cumsum, sort, sum, searchsorted, newaxis, arange, sqrt, log2, log, power, ceil, prod, zeros, ones, concatenate
from numpy.random import rand
from pylab import imread, imshow, plot, show
from weblogolib import LogoData, LogoOptions, LogoFormat, png_formatter, unambiguous_dna_alphabet
from itertools import repeat, chain, izip
from bisect import bisect_left
import functools
import time
from dreme import *

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
Y, list of strings. dataset of sequences
numsubs, number of subsequences
theta_motif, motif PWM matrix guess
theta_background_matrix, background PWM matrix guess
lambda_motif, motif frequency guess

Output:
theta_motif, motif PWM matrix
theta_background_matrix, background PWM matrix
lambda_motif, motif frequency
"""
def Online_EM(Y, numsubs, theta_motif, theta_background_matrix, lambda_motif):
    W = theta_motif.shape[0]#get the length of the motif
    #X = getSubsequences(Y,W)#this step may need to be removed for Online to save RAM
    #subsequences are grouped by sequences for normalization purposes
    #Is = [[sequenceToI(xij) for xij in xi] for xi in X]#list of indicator matrices, same dimensions as X
    s1_1 = lambda_motif#the expected number of occurrences of the motif
    s1_2 = theta_motif#the matrix holding the expected number of times a letter appears in each position, motif
    s2_2 = theta_background_matrix#the matrix holding the expected number of times a letter appears in each position, background
    n = 0#the counter
    nstart = 100000000#when to start averaging
    N = len(Y)#number of sequences
    #reserve some memory. this was when each sequence had only one subsequence
    """fractions = zeros(N)
    pwms = zeros((N,W,4))
    backgrounds = zeros((N,4))
    """
    #prepare lists because we don't know how many subsequences we have in total
    fractions = [None]*numsubs#pre-allocate space
    distances = [None]*numsubs#pre-allocate space
    pwms = list()
    backgrounds = list()
    expectations = list()
    fractions_sum = 0#should be deleted
    pwms_sum = zeros((W,4))#should be deleted
    backgrounds_sum = zeros((W,4))#should be deleted
    seqindpairs = [[(seqind,subseqind) for subseqind in range(len(Y[seqind]) - W)] for seqind in range(N)]
    seqindpairs = list(chain(*seqindpairs))
    random.shuffle(seqindpairs)
    #expectations.append(expected_LogLikelihood(Is, theta_motif, theta_background_matrix, lambda_motif))#add the expectation of the initial guess
    #p = Pool(20)
    truemotif = load('PWM_29_correct.npy')
    print "Running Online EM algorithm..."
    for seqindpair in seqindpairs:#iterate through each sequence index and start pair
        seqind = seqindpair[0]
        start = seqindpair[1]
        s = Y[seqind]#grab the whole sequence as a string
        L = len(s)#length of sequence
        #starts = range(0,L-W+1)
        #Is = [sequenceToI(s[k:k+W]) for k in starts]#indicator matrices of all subsequences in sequence s
        #random.shuffle(starts)#shuffle the starts for randomness
        #for start in starts:
            #I = Is[start]
        step = 0.05*pow(n+1,-0.6)#the online step size. For OLO6a
            #step = 0.025*pow(n+1,-0.6)#the online step size. For OLO6a
            #step = 1.0/10000
            #E-step
            #smooth
        left = max(0,start-W+1)
        right = min(L-W+1,start+W)
        middle = start - left
        Is = [sequenceToI(s[k:k+W]) for k in range(left,right)]
            #find expected Z values for all W-mers overlapping the current W-mer
        Z = [Z0_I(I,theta_motif, theta_background_matrix,lambda_motif) for I in Is]
            #find the Zs with parallel mapping
            #Z = p.map(functools.partial(Z0_I, theta_motif=theta_motif, theta_background_matrix=theta_background_matrix,lambda_motif=lambda_motif),Is[left:right])
        smooth(Z,W)#smooth the Z values
            #ds1_1 = Z0_I(I,theta_motif, theta_background_matrix,lambda_motif)
        ds1_1 = Z[middle]
            #print y
        I = Is[middle]
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
        fractions[n] = lambda_motif
        distances[n] = dist(theta_motif,truemotif)
            #pwms.append(theta_motif)
            #backgrounds.append(theta_background)
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
            #print n
            #the expected log likelihood, the objective function, based on current parameters
            #expectations.append(expected_LogLikelihood(Is, theta_motif, theta_background_matrix, lambda_motif))#add the expectation of the initial guess
    #print s1_2
    #x = load('NRF1_Motif.npy')
    #pylab.plot([dist(x,y) for y in pwms])
    #plot(fractions)
    #plot(distances)
    #show()
    #save('the_expectations', expectations)
    return theta_motif, theta_background_matrix, lambda_motif, fractions, distances
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
The smooth algorithm from Bailey and Elkan 1993. Constrains the sum of expectations in any window of width W to sum no
more than 1.

Input:
Z - a list of expected values. Its length is L - W + 1, where L is the length of the original sequence
W - the width of the motif

Output:
Z - the list of expected values normalized so that sums in any window of width W sum to 1 or less 
"""
def smooth(Z,W):
    for offset in range(W):
        for j in range(offset, len(Z) - W + 1, W):
            local = Z[j:j+W]
            localp = sum(local)
            if localp > 1:
                Z[j:j+W] = [l/localp for l in local]

"""
Similar to DREME's output_best_re function. Returns the best RE and number of sequence occurrences
in positive and negative set. 

Input:
re_pvalues - dictionary of p-values for all REs
pos_seqs - list of strings, positive sequence set
unerased_pos_seqs - list of strings, the original unerased positive sequence set
unerased_neg_seqs - list of strings, the original unerased negative sequence set
minw - minumum short motif width (default: 3)
maxw - maximum short motif width (default: 8)
ethresh - E-value threshold (default: 0.05)
log_add_pthresh - log10 of the p-value threshold in order to add to RE (default: log(0.01))
given_only - whether to search on negative strand (default: False)

Output:
best_re - RE object, the best RE found
pos - number of sequences in the positive set the best RE was found in
neg - number of sequences in the negative set the best RE was found in
best_log_pvalue - p-value of the best RE found
best_log_Evalue - E-value of the best RE found
unerased_log_Evalue - E-value of the best RE found using the unerased sequence sets
"""
def get_best_re(re_pvalues, pos_seqs, unerased_pos_seqs, unerased_neg_seqs, minw, maxw, ethresh, 
        log_add_pthresh, given_only):
    # get the best RE (lowest p-value) within width range
    candidates = [(re_pvalues[re][4], re) for re in re_pvalues if len(re)>=minw and len(re)<=maxw]
    if len(candidates) == 0:
        return("", "", 1e300, 1e300, 1e300)
    best_re = min(candidates)[1]
    r = re_pvalues[best_re]
    # used to allow 6 for consensus sequences but they shouldn't be created
    assert(len(r) == 5)
    pos = r[0]
    neg = r[2]
    best_log_pvalue = r[4]
    best_log_Evalue = best_log_pvalue + log(len(re_pvalues))
    # get the E-value if there had been no erasing
    unerased_log_pvalue = compute_exact_re_enrichment(best_re, unerased_pos_seqs, 
           unerased_neg_seqs, given_only)
    unerased_log_Evalue = unerased_log_pvalue[4] + log(len(re_pvalues))
    """
    # output the motif if significant
    if best_log_Evalue <= log(ethresh):
        pwm = make_pwm_from_re(best_re, pos_seqs, given_only=given_only)
        # disable timeout as now printing
        disable_timeout()
        # print the best RE
        write_xml_motif(xml_out, motif_num, best_re, pos, neg,
                best_log_pvalue, best_log_Evalue, unerased_log_Evalue, pwm, 
                get_matching_significant_components(best_re, re_pvalues, 
                        log_add_pthresh));
        # output a logo
        logo_out.output_logo(pwm, motif_num, best_re, False)
        # make rc motif
        pwm_rc = copy.deepcopy(pwm).reverseComplement()
        # output a logo
        logo_out.output_logo(pwm_rc, motif_num, get_rc(best_re), True)
    """
    return (best_re, pos, neg, best_log_pvalue, best_log_Evalue, unerased_log_Evalue)

"""
Similar to DREME's find_print function. Searches the sequences for the most significant RE.
Uses the default parameters from DREME.

Input:


pos_seqs - list of strings, positive sequence set
neg_seqs - list of strings, negative sequence set
unerased_pos_seqs - list of strings, the original unerased positive sequence set
unerased_neg_seqs - list of strings, the original unerased negative sequence set
ngen - beam width for generalization (default: 100)
nref - beam width for refinement (default: 1)
minw - minimum short motif width (default: 3)
maxw - maximum short motif width (default: 8)
mink - minimum width of core (default: 3)
maxk - maximum width of core (default: 8)
ethresh - E-value threshold (default: 0.05)
log_add_pthresh - log10 of the p-value threshold in order to add to RE (default: log(0.01))
given_only - whether to search on negative strand (default: False)
use_consensus - whether to convert REs longer than <maxk> to consensus sequence and refine (default: False)

Output:
best_word - RE object, the best RE found
pos - number of sequences in the positive set the best RE was found in
neg - number of sequences in the negative set the best RE was found in
best_log_pvalue - p-value of the best RE found
best_log_Evalue - E-value of the best RE found
unerased_log_Evalue - E-value of the best RE found using the unerased sequence sets
"""
def find_dreme_core(pos_seqs, neg_seqs, unerased_pos_seqs, unerased_neg_seqs, ngen=100, nref=1, minw=3, maxw=8, 
        mink=3, maxk=8, log_add_pthresh=log(0.01), ethresh=0.05, use_consensus=False, given_only=False):
    """
    Find a motif, print it, erase it.
    """
    #
    # Find core REs.
    #
    re_pvalues = re_find_cores(pos_seqs, neg_seqs, ngen, mink, maxk, 
            log_add_pthresh, given_only)

    #
    # Extend core REs to maximum width by finding new cores in flanking regions.
    #
    if (maxw > maxk):
        re_extend_cores(re_pvalues, ngen, mink, maxk, maxw, log_add_pthresh, 
                nref, use_consensus, pos_seqs, neg_seqs, given_only)

    #
    # Print the best word
    #
    (best_word, pos, neg, best_pvalue, best_Evalue, unerased_log_Evalue) = \
            get_best_re(re_pvalues, pos_seqs, unerased_pos_seqs, unerased_neg_seqs, minw, maxw, 
                ethresh, log_add_pthresh, given_only)
    _pv_format = "%3.1fe%+04.0f"
    rc_best_word = get_rc(best_word)
    pv_string = sprint_logx(best_pvalue, 1, _pv_format)
    ev_string = sprint_logx(best_Evalue, 1, _pv_format)
    unerased_ev_string = sprint_logx(unerased_log_Evalue, 1, _pv_format)
    print(("Best RE was {0:s} {1:s} p-value= {2:s} E-value= {3:s} "
               "Unerased_E-value= {4:s}").format(best_word, rc_best_word, 
                        pv_string, ev_string, unerased_ev_string))

    return (best_word, pos, neg, best_pvalue, best_Evalue, unerased_log_Evalue)

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
Generates a starting point by finding a sequence containing the core RE and follows the
heuristic in Bailey and Elkan to generate a starting point from the subsequence centered
on the instance of the core RE. Searches both strands.

Input:
core_re - RE of the motif core to search for in the sequences
pos_seqs - list of strings, the positive set of sequences
W - width of the pwm to generate
given_only - whether to search on both strands (default: False)

Output:
theta_motif - numpy array, the pwm generated from a subsequence containing the core motif
"""
def generate_starting_point(core_re, pos_seqs, W, given_only=False):
    ms = make_dna_re(core_re, given_only)#make a regular expression for both strands
    m = bisection(0.4,0.25,1)
    #search through each sequence until find a match far enough from the edges
    #search randomly
    inds = range(len(pos_seqs))
    random.shuffle(inds)
    for ind in inds:
        s = pos_seqs[ind]
        a = ms.search(s)#search sequence for RE
        if a:#if the sequence contains the RE
            L = len(s)
            (start, end) = a.span()
            core_width = end - start#width of the core found
            padding = (W - core_width)/2#padding to add to both sides (rounded down)
            start = start - padding
            end = end + padding
            if start >= 0 and end <= L:#check if starting point is actually inside sequence
                return startingPoint(sequenceToI(s[start:end]), m)

"""
The main online MEME algorithm.

Input:
Y, list of strings. dataset of sequences
W, width of motifs to search for
NPASSES, number of distinct motifs to search for

Output:
fractions
"""
def meme(Y,W,NPASSES):
    #6/28/13, check with initial conditions matching solution
    #p = Pool(64)
    #s=p.map(functools.partial(f,y=Y),range(64))
    #print s
    n = sum([max(0,len(y) - W + 1) for y in Y])#gets number of subsequences
    pos_seqs = Y
    print "Shuffling positive sequences..."
    neg_seqs = [ shuffle.dinuclShuffle(s) for s in pos_seqs ]#generate the negative sequence set by shuffling
    unerased_pos_seqs = copy.deepcopy(pos_seqs)
    unerased_neg_seqs = copy.deepcopy(neg_seqs)
    _dna_alphabet = 'ACGT'
    dprobs = get_probs(neg_seqs, _dna_alphabet)
    #generate first order Markov background based on nucleotide frequencies
    theta_background = array([[dprobs['A'], dprobs['C'], dprobs['G'], dprobs['T']]])
    lambda_motifs = list()
    theta_motifs = list()
    theta_background_matrices = list()
    fractionses = list()
    distanceses = list()
    for npass in range(NPASSES):
        #Use DREME to find a core seed
        (best_word, pos, neg, best_log_pvalue, best_log_Evalue, unerased_log_Evalue) = \
                find_dreme_core(pos_seqs, neg_seqs, unerased_pos_seqs, unerased_neg_seqs)
        #From DREME's best re, predict the fraction 
        lambda_motif = 1.0*pos/n/2#dividing by 2 accounts for the fact that RC not accounted for
        #From the best RE, search the positive sequence set and generate a starting PWM
        print 'Generating starting points from subsequences'
        theta_motif = generate_starting_point(best_word, pos_seqs, W)
        #theta_motif = load('NRSF_test.npy')
        theta_background_matrix = theta_background.repeat(theta_motif.shape[0],axis=0)#the initial guess for background is uniform distribution
        theta_motif, theta_background_matrix, lambda_motif, fractions, distances = Online_EM(Y, n, theta_motif, theta_background_matrix, lambda_motif)
        lambda_motifs.append(lambda_motif)
        theta_motifs.append(theta_motif)
        theta_background_matrices.append(theta_background_matrix)
        fractionses.append(fractions)
        distanceses.append(distances)
    return lambda_motifs, theta_motifs, theta_background_matrices, fractionses, distanceses
    
"""
Outputs the motif as a web logo. Saves fractions as .npy.

Input:
lambda_motif - a double, fraction of subsequences 
theta_motif - a numpy array, the PWM
theta_background_matrix - a numpy array, the background model
fractions - list of fractions at each step
k - the motif number index
outstr - the prefix for the output files
"""
def outputMotif(lambda_motif, theta_motif, theta_background_matrix, fractions, distances, k, outstr):
    _pv_format = "%3.1fe%+04.0f"
    f_string = sprint_logx(log(lambda_motif), 1, _pv_format)
    print(("Motif {0:s} had a fraction of {1:s}").format(str(k), f_string))
    print 'Saving motif as a png...'
    data = LogoData.from_counts(counts=theta_motif,alphabet=unambiguous_dna_alphabet)
    options = LogoOptions()
    options.title = 'Motif'
    forma = LogoFormat(data, options)
    fout = open(outstr + "Motif_" + str(k) + '.png', 'w')
    png_formatter(data, forma, fout)
    fout.close()
    print 'Saving fractions list'
    save(outstr + "Fractions_" + str(k), fractions)
    print 'Saving distances list'
    save(outstr + "Distances_" + str(k), distances)
    print 'Saving PWM'
    save(outstr + "PWM_" + str(k), theta_motif)
    #img = imread('results.png')
    #imshow(img)
    #show()
    #for now, just print, but will have to output a png later
    #print lambda_motif
    #print theta_motif
    #print theta_background_matrix
    
    

"""
The main executable function
"""
def main():
    usage = "usage: %prog [options] <input FASTA>"
    description = "The program applies the Online MEME algorithm to find motifs in a FASTA file"
    parser = ArgumentParser(description=description)
    parser.add_argument('fastafile', metavar='f', help='FASTA file containing the sequences')
    parser.add_argument("-p", "--processes", help="optional number of parallelized processes", type=int, default=1)
    parser.add_argument("-w", "--width", dest="width", help="Width of the motif to search for", type=int, default=0)
    parser.add_argument("-minw", dest="minwidth", help="Minimum width of the motif to search for", type=int, default=8)
    parser.add_argument("-maxw", dest="maxwidth", help="Maximum width of the motif to search for", type=int, default=40)
    parser.add_argument("-m", "--nummotifs", dest="nummotifs", help="Number of motifs to search for", type=int, default=1)
    parser.add_argument("-o", "--output", dest="output", help="Folder for all output files", default="output")    
    args = parser.parse_args()
    w = args.width
    if w == 0:
        minw = args.minwidth
        maxw = args.maxiwdth
    else:
        minw = w
        maxw = w
    nmotifs = args.nummotifs
    # make the directory (recursively)
    outdir = args.output
    outpre = outdir + "/"
    clobber = True
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not clobber:
                print >> sys.stderr, ("output directory (%s) already exists "
                "but DREME was not told to clobber it") % (outdir); sys.exit(1)
        else: raise
    print "Started at:"
    print time.ctime()
    #Use DREME's SeqIO to read in FASTA to list
    seqs = sequence.convert_ambigs(sequence.readFASTA(args.fastafile, None, True))
    lambda_motifs, theta_motifs, theta_background_matrices, fractionses, distanceses = meme(seqs,w,nmotifs)
    k = 1
    for lambda_motif, theta_motif, theta_background_matrix, fractions, distances in zip(lambda_motifs, theta_motifs, theta_background_matrices, fractionses, distanceses):
        outputMotif(lambda_motif, theta_motif, theta_background_matrix, fractions, distances, k, outpre)
        k = k+1
    print "Ended at:"
    print time.ctime()


if __name__=='__main__':
    main()