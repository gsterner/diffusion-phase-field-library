'''
Created on 28 dec 2010

@author: gustaf
'''
from ctypes import *
import calllib
import math
import random

GRID_SIZE = 256
XMAX      = 64.0
XMIN      = -64.0
DT        = 0.1    
TMAX_D      = range(2,22,2)
#TMAX_D      = [1]

TMAX_C       = [0.25, 1 , 10]
XMAX_C     = 128.0
XMIN_C     = -128.0
DT_C       = 0.05

TMAX_PFC   = 50
XMIN_P     = -102.4
XMAX_P     = 102.4
DT_P       = 0.01  
#EPSILON    = [1.0]  
EPSILON    = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

#---------------------------------------------------------------------------------

def start_diff(x):
    if (x < -0.5) or (x > 0.49999999):
        return 0;
    else:
        return 1;
#---------------------------------------------------------------------------------

def start_ch(x):
    if (x < 0):
        return 1;
    else:
        return -1;

#---------------------------------------------------------------------------------

def start_pfc(x):
    rnr = random.randint(0,10)
    return (rnr - 5) * 0.01

#---------------------------------------------------------------------------------

def diffusion(gs, xmx, xmn, dt, tmx):
    DoubleArray = c_double * GRID_SIZE
    gridSize = c_int(gs)
    max = c_double(xmx)
    min = c_double(xmn)
    dt = c_double(dt)
    tmax = c_double(tmx)
    out = DoubleArray()
    
    diffusion_c = calllib.get_diff_func(DoubleArray)

    diffusion_c(calllib.ALLOCATOR(start_diff), gridSize, max, min, dt, tmax, out)
    
    result = []
    for i in range(GRID_SIZE):
        result.append(out[i])
    return result
#---------------------------------------------------------------------------------

def cahnhilliard(gs, xmx, xmn, dt, tmx):
    DoubleArray = c_double * GRID_SIZE
    gridSize = c_int(gs)
    max = c_double(xmx)
    min = c_double(xmn)
    dt = c_double(dt)
    tmax = c_double(tmx)
    out = DoubleArray()
    
    cahnhilliard_c = calllib.get_ch_func(DoubleArray)

    cahnhilliard_c(calllib.ALLOCATOR(start_ch), gridSize, max, min, dt, tmax, out)
    
    result = []
    for i in range(GRID_SIZE):
        result.append(out[i])
    return result

#---------------------------------------------------------------------------------

def phasefieldcrystal(gs, xmx, xmn, dt, tmx, eps):
    DoubleArray = c_double * GRID_SIZE
    gridSize = c_int(gs)
    max = c_double(xmx)
    min = c_double(xmn)
    dt = c_double(dt)
    tmax = c_double(tmx)
    e = c_double(eps)
    out = DoubleArray()
    
    phasefieldcrystal_c = calllib.get_pfc_func(DoubleArray)

    phasefieldcrystal_c(calllib.ALLOCATOR(start_pfc), gridSize, max, min, dt, tmax, e, out)
    
    result = []
    for i in range(GRID_SIZE):
        result.append(out[i])
    return result


#---------------------------------------------------------------------------------

def xvals(gs, xmx, xmn, start_func):
    DoubleArray = c_double * GRID_SIZE
    gridSize = c_int(gs)
    max = c_double(xmx)
    min = c_double(xmn)
    out = DoubleArray()
    
    xvals_c = calllib.get_x_func(DoubleArray)

    xvals_c(calllib.ALLOCATOR(start_func), gridSize, max, min, out)
    
    result = []
    for i in range(GRID_SIZE):
        result.append(out[i])
    return result


#---------------------------------------------------------------------------------

def run():
    rArray = []
    for t in TMAX_D:
        rArray.append( diffusion(GRID_SIZE, XMAX, XMIN, DT, t))
    return rArray

#---------------------------------------------------------------------------------

def run_ch():
    rArray = []
    for t in TMAX_C:
        rArray.append( cahnhilliard(GRID_SIZE, XMAX_C, XMIN_C, DT_C, t))
    return rArray    

#---------------------------------------------------------------------------------
def run_pfc():
    rArray = []
    for e in EPSILON:
        rArray.append( phasefieldcrystal(GRID_SIZE, XMAX_P, XMIN_P, DT_P, TMAX_PFC, e))
    return rArray  

#---------------------------------------------------------------------------------

def getxvals():
    return xvals(GRID_SIZE, XMAX, XMIN, start_diff)

#---------------------------------------------------------------------------------

def getxvals_ch():
    return xvals(GRID_SIZE, XMAX_C, XMIN_C, start_ch)
#---------------------------------------------------------------------------------

def getxvals_pfc():
    return xvals(GRID_SIZE, XMAX_P, XMIN_P, start_pfc)

#---------------------------------------------------------------------------------
def amplitudes(M):
    amp = []
    for vec in M:
        amp.append(max(vec) - min(vec))
    return amp

#---------------------------------------------------------------------------------
def ch_analyt(xarray):
    retval = []
    onosqrttwo = 1/math.sqrt(2)
    for x in xarray:
        val = -1* math.tanh(x * onosqrttwo)
        retval.append(val)
    return retval

#---------------------------------------------------------------------------------
def diff_analyt(xarray, t):
    retval = []
    onofourt = 1/(4.0*t)
    onofpisq = 1/math.sqrt(4.0 * math.pi * t)
    for x in xarray:
        val = math.exp(-1*x*x*onofourt) * onofpisq
        retval.append(val)
    return retval
#---------------------------------------------------------------------------------
def pfc_analyt(epsArray):
    retval = []
    for ep in epsArray:
        val = 4*ep/3
        retval.append(math.sqrt(val))
    return retval
    
