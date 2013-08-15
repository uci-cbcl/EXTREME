OnlineMEME
==========

An Online implementation of the MEME algorithm

-----
Notes
-----
8/14/13:
Finally, after much crying, I am able to call the calc_entropy function with Cython,
as confirmed when it printed Entropy! when I used the calc_ent function in the MEME
class I made. After messing around with pointers and stuff, I can get C code to print
what I did in python and I can get python to print what I did in C code.

8/13/13:
Slowly, but surely, adding Cython support. Trying to wrap the likelihoods.c function so that I don't have to bother with the details
and I can just get an E-value. EXTERN and DEXTERN definitions are messing things up for some reason.

8/12/13:
Begin interfacing with C. E-value functions are in likelihood.c (about line 139).

8/9/13:
May want to consider trying without smoothing. Tested a lot. Motifs seem to coalesce
into a finite number of motifs. May consider just erasing if this is the case.

8/8/13:
Pregenerates all indicator matrices now. Next time, keep gamma_0=0.05 and alpha=0.6,
and try changing bisection gamma to 0.1 and see if I can get the same results (6/10 successes)
as bisection gamma=0.4.

8/7/13:
Shuffling subsequences instead of sequences now. Need to optimize it later. For now,
run multiple times and see how well it works.

8/6/13:
Added multiple outputs to a folder option. Save distances to a PWM for plotting later.
On the GM12878 set, it succeeds 6/10 times. It seems to first dip and then rise. 
The reason for failures may be that the algorithm has not seen enough instances
of the motif and so the fraction just keeps dipping until it hits a point of no
return.

8/2/13:
Removed plotting and added print statements so that I can run nohup.

8/1/13:
Added a DREME seeder. It took a few tries, but at a motif width of 20 managed to capture
the correct motif. Must find a way to terminate the program early in case of failed run.

7/27/13:
Extended the NRSF test by 10 base pairs on the left side. As expected, the PWM was still captured by the right side of the extended
PWM. Implemented DREME's method of storing sequences and getting nucleotide frequencies, no longer need pygr. Plotted distances.
For some reason, I captured a different motif from the sorted top 600 set. This motif had a higher fraction than the canonical
motif. Maybe because I used a more accurate background model. Stored in distances pngs. Idea for combining with DREME:
In each step, DREME finds a core, seed from instances of core using batch MEME's heuristic, use the DREME occurrences for the motif
fraction guess, Online MEME searches for maximum size motif, shrinks down or searches again for smaller
motifs, and then erase motif occurrences using Bayesian decision theory. Repeat this step until desired number of motifs
found or E-value threshold reached. 

7/26/13:
Hold off from speeding up DREME for now. Just focus on doing seed heuristic for now, because if seed heuristic doesn't work
then I just wasted my time. Extend the NRSF guess next time and see if that affects the middle still.

7/18/13:
Ran it on the intersected dataset with gamma_0 = 0.05. ~5,000,000 subsequences. Takes over an hour to complete. Got the correct
motif. I got this for the expected counts (I should do some math to figure out the best empirical gamma_0):

[[  3.66816003e-05   1.94069486e-05   5.85885736e-05   7.33175254e-05]
 [  1.73989938e-05   3.01240380e-05   2.62699048e-05   1.13888180e-04]
 [  1.88402226e-05   1.54034106e-04   1.21425579e-06   1.32706076e-05]
 [  1.68591485e-04   1.61766759e-10   1.74214575e-05   3.80694610e-06]
 [  2.28723451e-06   3.18960952e-06   1.76878963e-04   6.95718785e-06]
 [  8.95199351e-06   1.06088978e-04   5.63539555e-05   1.74481147e-05]
 [  1.88514567e-04   5.16481097e-17   3.19589787e-11   5.33415840e-14]
 [  9.03307774e-20   1.87958160e-04   5.32343674e-10   3.31747027e-18]
 [  1.86073636e-06   1.72447466e-04   8.24296119e-11   1.47330816e-05]
 [  8.45917941e-05   2.54954558e-05   2.61629141e-05   5.43548999e-05]
 [  3.50229811e-05   3.58553275e-05   2.09632097e-05   9.87766311e-05]
 [  8.26785637e-06   1.62108699e-23   1.82425288e-04   2.39415533e-21]
 [  3.91008458e-06   5.53666696e-13   1.86857366e-04   3.40091724e-21]
 [  1.54144648e-04   1.23110746e-05   1.89139208e-05   6.05309627e-06]
 [  1.36011811e-06   1.47828978e-04   3.84186626e-05   3.92587252e-06]
 [  1.92103503e-04   4.48274791e-16   2.67576344e-22   4.01529805e-20]
 [  1.99930234e-06   3.40984542e-07   1.89538657e-04   3.80471032e-18]
 [  2.79240466e-05   1.17481922e-04   2.53292311e-05   2.13611953e-05]
 [  4.63760590e-05   8.58640325e-06   8.34644507e-05   5.39810773e-05]
 [  2.79646290e-05   1.03917789e-04   3.56657171e-05   2.51252625e-05]
 [  2.60511430e-05   1.16332327e-04   1.04380186e-05   4.00787910e-05]]

lambda = 0.000242620408426

7/17/13:
Implemented smoothing for Online. It calculates the expected Z value for all overlapping W-mers and then performs the smoothing
algorithm on this subsequence. Tested this new feature on the top 600 NRSF K562 set and was able to achieve the correct motif even
when gamma_0 = 0.1. However, it still failed at gamma_0 = 0.2. With smoothing, the program takes about 8 minutes to complete
the 600 sequence dataset (~500,000 subsequences).

7/16/13:
Tried online MEME on the Myers data. Only looked at rep 1 peaks that overlapped with rep 2 peaks for the K562 NRSF dataset. For the top
600 peaks, without any trimming, Online MEME usually converges to the correct motif. If gamma_0 is too high, it will converge to 
repetitive elements or converge to "pitfalls". If gamma_0, it will not converge fast enough and the pwm will not be very "crisp".
This is even more of a problem when I included every single overlapping peak. For top 600, gamma_0 of 0.05 works fine. It has
about 500,000 subsequences and takes 45 seconds to run. For the entire overlapping dataset, gamma_0 of 0.01 works okay, but 
clearly does not converge fast enough. If gamma_0 gets too big, it will converge to a repetitive motif. I could fix this by either 
implementing the smoothing feature or running through the dataset multiple times. 

7/15/13:
MEME-chip sucks. Downloaded K562 NRSF (1st replicate) broadpeaks file from Myers lab ENCODE as a real world dataset. gamma_0=0.2
seems to be the sweet spot. Used a simulated dataset of NRSF, mean length 100, 80% of sequences contain motif, 1000 sequences.
Online MEME worked very fast, but plots appear "jagged", even when subsequences were shuffled. I have to figure out how to do
Online smoothing next.


7/13/13:
Added the smooth function to the Batch MEME algorithm. Also added the weblogo display to Batch MEME. Batch MEME displays the
plot of expected log likelihoods at the end of the Batch MEME algorithm. Tested a dataset with 100 sequences, 80% of which contain
the NRSF motif. The official MEME and the Batch MEME algorithm correctly found the motif. Must test this on Online MEME next.
Then will add smoothing feature to Online MEME and account for reverse complements.

7/8/13:
The algorithms can now save and plot the expected log likelihood at each iteration. As expected, both algorithms
asymptotically reach a maximum, although only the batch algorithm reaches the global maximum. Tested on a NRF1
mock data set, 5000 sequences of length 10, 50% motif fraction. Also corrected some spelling errors.

7/5/13:
Added long sequence functionality to Online-MEME. Get same results as batch mode for 12-bp long NRF1 example. However,
Batch MEME, Online MEME, and even the web-based MEME fails for longer NRF1 sequences. Also changed the way parameters are stored.
Next step is to plot the expected-likelihood and account for reverse complements.

7/2/13:
Tested the efficacy of averaging. Confirmed what the Online EM paper said. Smaller initial step removes bias, lowers variance,
but slows convergence. Averaging can lower the variance for bigger initial steps, but requires memory to store past parameters.
General rule of thumb: only use averaging if you do not have that many samples. If you have a lot of samples, then just use 
a small initial step, since the Online EM algorithm will converge to the true parameters before it reaches the last observation.
Next step: modify Online MEME to handle multiple subsequences per sequence.

7/1/13:
The code works better if the data is actually randomized... May just use the OLO6 method to save memory. Seems to be a tradeoff
between convergence rate and accuracy/consistency.

6/29/13:
Changed the step size forumla. Works pretty well so far. Perturbed background model and it performed fairly well. Needs more testing.
It appears that OLO6a has superior performance to OLO6. For the current perturbed initial conditions, OLO6a converged better to correct
solution. Should do more parameter testing just to make sure.

6/28/13:
Finished writing the prototype for the Online MEME algorithm. Assumes the FASTA file holds sequences with one subsequence each.
Will test on the NRF1_TestMotif.fasta file.

6/27/13:
Tested the Batch MEME algorithm on a test set of data. The test data consisted of 10000 sequences. Each sequence
was 10 bp long. Half of the sequenced were generated from a uniform background. The other half were generated
from the NRF1 motif. The batch algorithm correctly called the motif and the fraction using the heuristic for
maximum expectation. The algorithm can not yet handle erasure or priors.
Added minimum width and number of motifs to find option to Batch MEME.

5/27/13:
Got the EM and heuristics part working. Need to figure out how to speed up
the heuristics part, way too slow. Perhaps use multiprocessing? Exported
the NRSF PWM as a numpy array. Should use this as a starting point. The current
starting points are not working well.

5/25/13:
Finished bulk of the main MEME algorithm. Must now work on output and actually
running it.

5/16/13:
Completed the random sampling of subsequences and relative entropy per column. Work on deriving starting points next.

5/14/13:
Implemented E and M steps. Updated functions to use onle the indicator matrices.Next time need to work on the random sampling part.

5/7/13:
timeit.timeit('px_background(s,theta_background)','from __main__ import px_background,s,theta_background', number=10000) took 0.2 s.
timeit.timeit('px_background(str(db[db.keys()[50]][50:60]),theta_background)','from __main__ import px_background,db,theta_background', number=10000) took 2 s.
timeit.timeit('x=getSubsequences(db,10)','from __main__ import db,getSubsequences', number=1) takes 1.13 s.
timeit.timeit('[[px_background(a,theta_background) for a in A] for A in x]','from __main__ import x,px_background,db,theta_background,getSubsequences', number=1) takes 0.38 s.
Clearly, pygr indexing is much slower than direct memory access. However, converting to strings is time intensive as well. Something to keep in mind for speed and memory purposes.
Based on my estimates, it will take twice as long reading from pygr each time compared to converting all to strings first. For the regular MEME version, I will load all into
memory, but I should keep in mind that for the Online version memory will be a constraint.
M00256.pfm is a motif file in the MotifList format corresponding to NRSF
M00652.pfm is a motif file in the MotifList format corresponding to NRF1
