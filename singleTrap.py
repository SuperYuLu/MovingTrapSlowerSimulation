# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:53:06 2017

@author: Yu Lu
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
        self._trapCenter = trapCenter
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
        frontCoilCenter = self.coilSpace / 2 + self._trapCenter
        backCoilCenter = -self.coilSpace / 2 + self._trapCenter
        B = 0
        for centerPos, layers, turns in zip( # for front and back coil 
                    [frontCoilCenter, backCoilCenter], 
                    [self.numLayersFront, self.numLayersBack], 
                    [self.numTurnsPerLayerFront, self.numTurnsPerLayerBack]):
            
            for l in range(layers):
                layerRadius = self.coilRadius + (l + 0.5) * self.wireDia # radius for specific layre 
                for t in range(turns):
                    windingCenter = centerPos  / 2 - (layers / 2  + 0.5) * self.wireDia + t * self.wireDia # center position for single winding
                    current = curr if centerPos > self._trapCenter else -curr
                    B += u0 * current * layerRadius**2 / 2 / (((z - windingCenter)**2 + layerRadius**2)**1.5)
        return B

        
    
    def plotField1D(self, pos, B):
        """
        Plot the magnetic field distribution along 
        a single dimention
        """
        try:
            assert len(pos) == len(B)
        except AssertionError:
            print("position and field has different length: ({:2d}, {:2d})".format(len(pos), len(B)))

        else:
            fig, ax = plt.subplots()
            ax.plot(pos*1e3, B, '*')
            plt.show()
        


class slower():
    def __init__(self, initialV, finalV, accRatio, current):
        """
        Initialize slower geometry and dynamics settings;
        All geometry parameters are calculating from trap center;
        """
        # Geometry Settings
        self.coilSpace = 10e-3
        self.trapSpace = self.coilSpace / 2.
        self.numTraps = 480
        self.divTrapNum = 180
        
        self.totalLength = (self.numTraps - 1) * self.trapSpace
        self.stage1Length = self.divTrapNum * self.trapSpace
        self.stage2Length = (self.numTraps - self.divTrapNum) * self.trapSpace

        # Dynamics Settings 
        self.initialV = initialV
        self.finalV = finalV
        self.accRatio = accRatio
        
        self.stage1Acc, self.stage2Acc = self.__accValue__()
        self.middleV = np.sqrt(finalV**2 - 2 * self.stage1Acc * self.stage1Length)
        
        self.stage1Time, self.stage2Time, self.totalTime = self.__time__()

        # Electronics 
        self.current = current
        
    def __accValue__(self):
        """
        Calculate the acceleration value for 1st and 2nd stage
        """
        
        a1 = (self.finalV**2 - self.initialV**2) / (2 * (self.stage1Length + self.accRatio * self.stage2Length))
        a2 = a1 * self.accRatio
        return a1, a2

    def __time__(self):
        """
        Calculate total time for deceleration
        """
        t1 = (self.middleV - self.initialV) / self.stage1Acc
        t2 = (self.finalV - self.middleV) / self.stage2Acc
        t = t1 + t2
        return t1, t2, t
    
    def __str__(self):

        print("\nSlower Parameters:\n================================")
        print("Num Traps:  {:<25} \t Length[m]:  {:<}".format(self.numTraps, self.totalLength),\
              "\nInitial V[m/s]:  {:<25} \t Final V[m/s]:  {:<}".format(self.initialV, self.finalV),\
              "\nStage1 Acc[m/s2]:  {:<25.2f} \t Stage2 Acc[m/s2]:  {:<.2f}".format(self.stage1Acc, self.stage2Acc),\
              "\nStage1 time[ms]:  {:<25.2f} \t Stage2 time[ms]:  {:<.2f}".format(self.stage1Time*1e3, self.stage2Time*1e3),\
              )
        return '0'
        
    def effectiveOnAxisMagField(self, curr, z):
        """
        Calculate the effecitive magnetic field on coil axis in the co-moving
        frame of the trap it self
        """
        return 0
            
        
    
