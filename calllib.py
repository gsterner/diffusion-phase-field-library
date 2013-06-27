'''
Created on 28 dec 2010

@author: gustaf
'''

from ctypes import *

ALLOCATOR = CFUNCTYPE(c_double, c_double)

def get_diff_func(DoubleArray):
    difflib = cdll.diffl
    diffunc = difflib.diffusion
    diffunc.restype = c_int
    diffunc.argtypes  = [ALLOCATOR, c_int, c_double, c_double, c_double, c_double, DoubleArray]
    return diffunc
    
def get_x_func(DoubleArray):
    difflib = cdll.diffl
    xvaluesfunc = difflib.xvalues
    xvaluesfunc.restype = c_int
    xvaluesfunc.argtypes  = [ALLOCATOR, c_int, c_double, c_double, DoubleArray]
    return xvaluesfunc

def get_ch_func(DoubleArray):
    difflib = cdll.diffl
    chfunc = difflib.cahnhilliard
    chfunc.restype = c_int
    chfunc.argtypes  = [ALLOCATOR, c_int, c_double, c_double, c_double, c_double, DoubleArray]
    return chfunc

def get_pfc_func(DoubleArray):
    difflib = cdll.diffl
    pfcfunc = difflib.phasefieldcrystal
    pfcfunc.restype = c_int
    pfcfunc.argtypes  = [ALLOCATOR, c_int, c_double, c_double, c_double, c_double, c_double, DoubleArray]
    return pfcfunc