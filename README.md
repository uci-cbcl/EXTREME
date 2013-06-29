OnlineMEME
==========

An Online implementation of the MEME algorithm

-----
Notes
-----
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
