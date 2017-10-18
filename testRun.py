# testRun.py --- 
# 
# Filename: testRun.py
# Description: 
# 
# Author:    Yu Lu
# Email:     yulu@utexas.edu
# Github:    https://github.com/SuperYuLu 
# 
# Created: Mon Oct 16 07:29:14 2017 (-0500)
# Version: 
# Last-Updated: Mon Oct 16 07:32:51 2017 (-0500)
#           By: yulu
#     Update #: 1
# 

from singleTrap import *
import numpy as np 

trapNum = np.array([x for x in range(1,480)])
initialV = 500
finalV = 0
accRatio = 1
current = 500

N = 10
p = singleTrap(N,
               initialV,
               finalV,
               accRatio,
               current)

z, B, B_eff = p.onAxisMagField()
p.plotField1D(z, B, B_eff)
print(p.trapFieldCenter(z, B))
print(p.fieldPeak(z, B_eff))
