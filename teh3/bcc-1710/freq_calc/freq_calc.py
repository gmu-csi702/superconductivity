#!/bin/env python3

import os, sys, path
import numpy as np

bv = np.loadtxt("bulkmodfit.out")
vol = (bv[:,1])**(1.0/6.0)
bm = (bv[:,2])**(1.0/2.0)

#calculate c constant
w = np.zeros((len(bm)))
w_h3cl56 = 1338
x_h3cl56 = 527
h_h3cl = 1578
c_w_h3cl = w_h3cl56* ( 1.0/(bm[0]*vol[0]) )
c_x_h3cl = x_h3cl56* ( 1.0/(bm[0]*vol[0]) )

table = []
i = 0
for i in range(len(bm)):
    w = c_w_h3cl * bm[i] * vol[i]
    s = c_x_h3cl * bm[i] * vol[i]
    
    x = np.array((w, h_h3cl, s))
    table.append(x)
    
np.savetxt("h3cl_freq", table, fmt='%.2f')
    

w_h3s56 =  1560
x_h3s56 = 615
h_h3s = 1840
c_w_h3s = w_h3s56* ( 1.0/(bm[0]*vol[0]) )
c_x_h3s = x_h3s56* ( 1.0/(bm[0]*vol[0]) )

table = []
i = 0
for i in range(len(bm)):
    w = c_w_h3s * bm[i] * vol[i]
    s = c_x_h3s * bm[i] * vol[i]
    
    x = np.array((w, h_h3s, s))
    table.append(x)
    
np.savetxt("h3s_freq", table, fmt='%.2f')
