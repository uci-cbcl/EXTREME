cdef extern from "mtype.h":
    cdef enum MOTYPE:
        Tcm, Oops, Zoops

cdef extern from "meme.h":
    ctypedef struct DATASET:
        int nmotifs
    ctypedef struct MODEL:
        int min_w

    void calc_entropy(MODEL *model, DATASET *dataset)

cdef class MEME:
    cdef MODEL* m
    cdef DATASET* d
    def __cinit__(self, theta_motif, theta_background, lambda_motif, sequences):
        cdef MODEL mo
        mo.min_w = 1
        self.m = &mo 
        x = 1

    def calc_ent(self):
        calc_entropy(self.m, self.d)
