# functions.py --- 
# 
# Filename: functions.py
# Description: 
#         A collections of functions for moving
#     trap slower calculation & simulation
# Author:    Yu Lu
# Email:     yulu@utexas.edu
# Github:    https://github.com/SuperYuLu 
# 
# Created: Wed Oct 18 21:56:47 2017 (-0500)
# Version: 
# Last-Updated: Wed Oct 18 22:38:32 2017 (-0500)
#           By: yulu
#     Update #: 9
# 


import numpy as np

u0 = 4 * np.pi * 1e-7
mj = 0.5
gj = 2
ub = 9.274e-24 # [J/T]
amu = 1.66e-27 # [Kg] Litium mass        
kb = 1.38e-23 
def trapDepth(B):
    """
    Calculate the trap depth of a magentic trap for 
    element with mass = [mass] amu, return the trap
    depth in terms of temperature in Kelvin
    """
    T  = 2 * ub * gj * mj * abs(B) / kb
    return T


