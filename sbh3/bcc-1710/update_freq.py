#!/bin/env python3

import sys,os
import numpy as np


def freq_calc(dos_file, bmv_file, w56, x56, h56, outfilename ):
    vals = np.loadtxt(dos_file)
    bv = np.loadtxt(bmv_file)
    vol = (bv[:,1])**(1.0/6.0)
    bm = (bv[:,2])**(1.0/2.0)

    #calculate c constant
    c_w = w_56 * ( 1.0/(bm[0]*vol[0]) )
    c_x = x_56 * ( 1.0/(bm[0]*vol[0]) )
    
    table = []
    for i in range(len(bm)):
        w = c_w * bm[i] * vol[i]
        x = c_x * bm[i] * vol[i]
        
        f = np.array((w, h56, x))
        table.append(x)


    np.savetxt(outfilenname, table, fmt='%.2f')



