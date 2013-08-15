#cdef extern from "logs.h":
#    void init_log()

#init_log()

print 'wassaup'

cdef extern from "macros.h":
    ctypedef int BOOLEAN 

cdef extern from "mtype.h":
    cdef enum MOTYPE:
        Tcm, Oops, Zoops

cdef extern from "meme.h":
    ctypedef struct SAMPLE:
        long length
    ctypedef struct DATASET:
        int alength
        int nmotifs
        char *alphabet
        double *back
    ctypedef struct MODEL:
        int w
        MOTYPE mtype
        BOOLEAN pal
        BOOLEAN invcomp
        double *rentropy
        double ic
        double logpv
        double logev
        double llr

    void calc_entropy(MODEL *model, DATASET *dataset)
    void dq_test(MODEL *model, DATASET *dataset)

cdef class MEME:
    cdef MODEL m
    cdef DATASET d
    def __cinit__(self, theta_motif, theta_background, lambda_motif, sequences):
        #cdef MODEL mo
        #cdef DATASET da
        self.m.w = 2
        self.m.mtype = Tcm
        self.m.pal = 0
        #self.m = &mo 
        self.d.alength = 4
        self.d.alphabet = "ACGT"
        #self.d = &da
        x = 2

    def calc_ent(self):
        calc_entropy(&self.m, &self.d)
        
    def dq(self):
        dq_test(&self.m,&self.d)
        print self.m.ic
        
