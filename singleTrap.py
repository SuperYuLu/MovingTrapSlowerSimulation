# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:53:06 2017

@author: Beam
"""
import numpy as np
import matplotlib.pyplot as plt 

# Constants needed 
u0 = 4 * np.pi * 1e-7

class singleTrap:
    def __init__(self, trapCenter = 0):
        """
        Initialize trap coil geometry
        """
        self.trapCenter = trapCenter
        self.coilRadius = 5.1e-3 # Measure to the inner layer of coil
        self.coilSpace = 10e-3 # Measure from the center of front/back coil
        self.wireDia = 0.405e-3 # AWG 26
        self.numLayersFront = 4
        self.numTurnsPerLayerFront = 4
        self.numLayersBack = 4
        self.numTurnsPerLayerBack = 2
    
    def onAxisMagField(self, curr, z):
        """
        Calculate the magnetic field generated from the anti-helmholtz coil
        under current "curr" at position z (from the trap center), in the 
        labtory frame
        """
        frontCoilCenter = self.coilSpace / 2 + self.trapCenter
        backCoilCenter = -self.coilSpace / 2 + self.trapCenter
        B = 0
        for centerPos, layers, turns in zip( # for front and back coil 
                    [frontCoilCenter, backCoilCenter], 
                    [self.numLayersFront, self.numLayersBack], 
                    [self.numTurnsPerLayerFront, self.numTurnsPerLayerBack]):
            
            for l in range(layers):
                layerRadius = self.coilRadius + (l + 0.5) * self.wireDia # radius for specific layre 
                for t in range(turns):
                    windingCenter = centerPos  / 2 - (layers / 2  + 0.5) * self.wireDia + t * self.wireDia # center position for single winding
                    current = curr if centerPos > self.trapCenter else -curr
                    B += u0 * current * layerRadius**2 / 2 / (((z - windingCenter)**2 + layerRadius**2)**1.5)
        return B
    
    
        
            
                    
                
            
        
    