#!/usr/bin/env python

import os, re, sys,math

from numpy import *

# frequency of envt 2
#freq_patch_2 = list(arange(0.01,0.99,0.02))
freq_patch_2 = list(arange(0.01,0.99,0.02))

# avarage switch rate
#sbar = list(arange(-1.5, 0.5, 2.0/50))
sbar = list(arange(-1.5, 0.5, 2.0/50))
d = [ 0.1, 0.5, 1.0 ]

exe = "./numsolve"

max_iter = 1000000

# initial values for patch freqs and relatedness coeffs
fval_init = " ".join([ str(1.0/8) for xi in range(0,8) ])

# initial values for the reproductive values
vval_init = " ".join([ str(1.0) for xi in range(0,8) ])

ctr = 0

for f2 in freq_patch_2:
    for sbar_i in sbar:
        for d_i in d:
            s2 = sqrt(((1.0-f2)/f2) * 10**(2*sbar_i))
            s1 = 10**(2*sbar_i) / s2
            print("echo " + str(ctr))
            print(exe + " " + str(max_iter) + " " 
                    + str(2) + " " 
                    + str(1) + " " 
                    + str(d_i) + " " 
                    + str(s1) + " " 
                    + str(s2) + " " 
                    + str(0.5) + " " 
                    + str(0.5) + " " 
                    + str(0.5) + " " 
                    + str(0.5) + " " 
                    + fval_init + " "
                    + vval_init + " ")
            ctr+=1
        
