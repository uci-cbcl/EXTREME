import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc
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
    ctypedef double **THETA
    ctypedef struct SAMPLE:
        long length
    ctypedef struct DATASET:
        int alength
        int n_samples
        char *alphabet
        double *back
        SAMPLE **samples
    ctypedef struct MODEL:
        THETA obs
        int w
        int nsites_dis
        MOTYPE mtype
        BOOLEAN pal
        BOOLEAN invcomp
        double *rentropy
        double rel
        double ic
        double logev
        double llr

    void calc_entropy(MODEL *model, DATASET *dataset)
    void dq_test(MODEL *model, DATASET *dataset)

cdef double** npy2d_double2d(np.ndarray[double, ndim=2, mode="c"] a):
    cdef double **tmp_c = <double**> malloc(a.shape[0] * sizeof(double*))
    for k in range(a.shape[0]):
        tmp_c[k] = &a[k,0]
    return tmp_c

cdef class MEME:
    cdef MODEL m
    cdef DATASET d
    cdef pwm
    cdef background
    def __cinit__(self, np.ndarray[double, ndim=2, mode="c"] theta_motif, np.ndarray[double, ndim=1, mode="c"] theta_background, lambda_motif, sequences, int nsites_dis):
        #these lines make sure the arrays don't disappear
        self.pwm = theta_motif
        self.background = theta_background
        self.m.w = theta_motif.shape[0]
        self.m.mtype = Tcm
        self.m.pal = 0
        self.m.invcomp = 1
        self.m.nsites_dis = nsites_dis
        #do not need to reserve space for rentropy. C code already did that
        #reserve some space for rentropy, which will get assigned in C
        #cdef double *rentropy = <double*> malloc(self.m.w * sizeof(double))
        #self.m.rentropy = rentropy #<double*> malloc(self.m.w * sizeof(double))
        #print theta_motif
        #print theta_background
        #malloc tricks to convert 2D numpy array to double**
        self.m.obs = npy2d_double2d(theta_motif)
        self.d.back = &theta_background[0]
        self.d.alength = len(theta_background)
        self.d.alphabet = "ACGT"
        self.d.n_samples = len(sequences)
        #do some malloc tricks to get an array of SAMPLE pointers
        cdef SAMPLE* ss = <SAMPLE*> malloc(self.d.n_samples * sizeof(SAMPLE))
        self.d.samples = <SAMPLE**> malloc(self.d.n_samples * sizeof(SAMPLE*))
        for i in range(self.d.n_samples):
            ss[i].length = len(sequences[i])
            self.d.samples[i] = &ss[i]
        x = 3

    def calc_ent(self):
        calc_entropy(&self.m, &self.d)
    
    def get_logev(self):
        return self.m.logev
    
    def get_rentropy(self):
        r = np.zeros(self.m.w)
        for i in range(self.m.w):
            r[i] = self.m.rentropy[i]
        return r
        
    def dq(self):
        dq_test(&self.m,&self.d)
        print self.m.ic
        print self.pwm
        print self.background
        

        
        
