OnlineMEME
==========

An Online implementation of the MEME algorithm

-----
Notes
-----
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
