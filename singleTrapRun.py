# singleTrapRun.py --- 
# 
# Filename: testRun.py
# Description: 
#            Run single trap analysis based on the
#          given parameters in settings.py
# Author:    Yu Lu
# Email:     yulu@utexas.edu
# Github:    https://github.com/SuperYuLu 
# 
# Created: Mon Oct 16 07:29:14 2017 (-0500)
# Version: 
# Last-Updated: Wed Oct 18 23:09:58 2017 (-0500)
#           By: yulu
#     Update #: 25
# 

from settings import *
from movingTraps import *
from functions import trapDepth
import numpy as np 

def main():
    global sigleTrapRunNum, initialV, finalV, accRatio, current

    p = singleTrap(singleTrapRunNum,
                   initialV,
                   finalV,
                   accRatio,
                   current)

    z, B, B_eff = p.onAxisMagField()
    
    print("\n===============================")
    print("Trap number: {:}".format(singleTrapRunNum))
    print("Trap field center at {:.2f} mm".format(p.trapFieldCenter(z, B) * 1e3))
    print("Trap field peak: front {:.2f}T  back {:.2f}T".format(*p.fieldPeak(z, B_eff)))
    print("Trap depth: {:.2f} mK".format(trapDepth(min(p.fieldPeak(z, B_eff))) * 1e3 ))
    p.plotField1D(z, B, B_eff)
    
if __name__ == "__main__":
    main()


