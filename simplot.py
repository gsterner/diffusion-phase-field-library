'''
Created on 28 dec 2010

@author: gustaf
'''

import pylab as PL
import simrun
import time

#===============================================================================

D = simrun.run()
x = simrun.getxvals()
PL.figure(1)
labels = []
i = 0
for d in D:
   times = simrun.TMAX_D
   a = simrun.diff_analyt(x, times[i])    
   PL.plot(x,d, 'o')
   PL.plot(x,a)
   i += 1
PL.legend()
PL.title("Diffusion Equation, different simulation times. Dots - Numerical. Lines - Analytical")
PL.xlabel('x')
PL.ylabel('PSI')
PL.show()
#===============================================================================

#------------------------------------------------------
#C = simrun.run_ch()
#xc = simrun.getxvals_ch()
#chanalyt = simrun.ch_analyt(xc)
#PL.figure(1)
#labels = []
#i = 0
#for c in C:
#    times = simrun.TMAX_C
#    ct = 't = ' + str(times[i])
#    labels.append(ct)
#    PL.plot(xc,c, '-o', label = ct)
#    i += 1
#PL.plot(xc,chanalyt, label = 'Analytical')
#PL.legend()
#PL.axis([-7, 7, -1.3, 1.3])
#PL.title("Cahn-Hilliard with different simulation times")
#PL.xlabel('x')
#PL.ylabel('PSI')
#PL.show()

#------------------------------------------------------

#===============================================================================
# tic = time.clock()
# P = simrun.run_pfc()
# toc = time.clock()

# print '>>> TIME Elapsed:', toc - tic

# xp = simrun.getxvals_pfc()
# PL.figure(2)
# PL.plot(xp,P[len(P)-1], 'o-')
# PL.show()

# amps = simrun.amplitudes(P)
# epss = simrun.EPSILON
# analamps = simrun.pfc_analyt(epss)

# PL.figure(1)
# PL.plot(epss,amps, 'o', label = 'Numerical')
# PL.plot(epss,analamps, label = 'Analytical')
# PL.legend()
# PL.title('Phase Field Crystal Equation solutions. Amplitude as a function of epsilon')
# PL.xlabel('espilon')
# PL.ylabel('Amplitude')
# PL.show()
#===============================================================================

# step = (simrun.XMAX_P - simrun.XMIN_P)/simrun.GRID_SIZE
# x_start = PL.arange(simrun.XMIN_P, simrun.XMAX_P,step )
# y_start = []
# for x in x_start:
    # y_start.append(simrun.start_pfc(x))
    
PL.plot(x_start, y_start)
# print PL.mean(y_start)
