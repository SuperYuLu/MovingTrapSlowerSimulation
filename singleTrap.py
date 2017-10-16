# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:53:06 2017

@author: Yu Lu
"""
import numpy as np
import matplotlib.pyplot as plt 

# Constants needed 
u0 = 4 * np.pi * 1e-7
mj = 0.5
gj = 2
ub = 9.274e-24 # [J/T]
M = 1.66e-27 * 7 # [Kg] Litium mass        


class slower:
    
    ###########
    # GLOBALS #
    ###########
    
    ## GEOMETRY OF A SINGLE TRAP COIL
    ################################
    coilRadius = 5.1e-3 # Measure to the inner layer of coil
    coilSpace = 10e-3 # Measure from the center of front/back coil
    wireDia = 0.405e-3 # AWG 26
    numLayersFront = 4
    numTurnsPerLayerFront = 4
    numLayersBack = 4
    numTurnsPerLayerBack = 2
    
    
    ## GEOMETRY OF SLOWER
    #####################
    trapSpace = coilSpace / 2.
    numTraps = 480
    divTrapNum = 180
    
    totalLength = (numTraps - 1) * trapSpace
    stage1Length = divTrapNum * trapSpace
    stage2Length = (numTraps - divTrapNum) * trapSpace
    middlePos = trapSpace * (divTrapNum - 1)
    
    def __init__(self, initialV, finalV, accRatio, current):
        """
        Initialize slower geometry and dynamics settings;
        All geometry parameters are calculating from trap center;
        """
        
        ## DYNAMICS OF SLOWER
        #####################
        
        self.initialV = initialV
        self.finalV = finalV
        self.accRatio = accRatio
        
        self.stage1Acc, self.stage2Acc = self.__accValue__()
        self.middleV = np.sqrt(finalV**2 - 2 * self.stage1Acc * self.stage1Length)
        self.stage1Time, self.stage2Time, self.totalTime = self.__time__()

        ## ELECTRONICS OF SLOWER
        ########################
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
              "\nTotal time[ms]:  {:<25.2f} \t Middle V[m/s]:  {:<.2f}".format(self.totalTime * 1e3, self.middleV)\
              )
        return '\n'

    def calcTrapVelocity(self, trapNum):
        """
        Calculate the effective moving velocity of specified trap
        """
        pos = self.trapSpace * (trapNum - 1)
        if trapNum <= self.divTrapNum:
            v = np.sqrt(self.initialV**2 + 2 * self.stage1Acc * pos)
        else:
            v = np.sqrt(self.middleV**2 + 2 * self.stage2Acc * (pos - self.middlePos))
        return v

    def calcTrapOnTime(self, trapNum):
        """
        Find trap turn on time and pulse length
        assumes trap at maximium depth when particle arrives at the center
        """
        
        if trapNum == 1:
            tPeriod = 2 * self.trapSpace / 2 / self.calcTrapVelocity(trapNum)
            tOn = -0.5 * tPeriod
        else:
            pos = self.trapSpace * (trapNum - 1)
            tOn_pre, tPeriod_pre = self.calcTrapOnTime(trapNum - 1)
            tOn = tOn_pre + 0.5 * tPeriod_pre
            tPeriod = 2 * self.trapSpace / 2 / self.calcTrapVelocity(trapNum)
        return tOn, tPeriod
            

    def calcTrapAcc(self, trapNum):
        """
        find  acceleration value for a given trap number
        """
        
        acc = self.stage1Acc if trapNum <= self.divTrapNum else self.stage2Acc
        return acc
        
    def effectiveOnAxisMagField(self, curr, z):
        """
        Calculate the effecitive magnetic field on coil axis in the co-moving
        frame of the trap it self
        """
        return 0
            
        
class singleTrap(slower):
    def __init__(self, trapNum, initialV = 480, finalV = 50, accRatio = 1, current = 400):
        """
        Initialize trap coil geometry
        """
        super().__init__(initialV, finalV, accRatio, current)
        self.trapCenter = self.trapSpace * (trapNum - 1)
        self.trapNum = trapNum
        
        self.trapVelocity = self.calcTrapVelocity(trapNum)
        self.trapTurnOnTime, self.pulseLength  = self.calcTrapOnTime(trapNum)
        self.trapAcc = self.calcTrapAcc(trapNum)
        
        
    def onAxisMagField(self, z = []):
        """
        Calculate the magnetic field generated from the anti-helmholtz coil
        under current "curr" at position z (from the trap center), in the 
        labtory frame
        """
        if z ==[]:
            z = np.linspace(-self.coilSpace*2, self.coilSpace * 2, 100) + self.trapCenter
        else:
            pass
        
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
                    windingCenter = centerPos - (turns / 2  + 0.5) * self.wireDia + t * self.wireDia # center position for single winding
                    current = self.current if centerPos > self.trapCenter else -self.current
                    B += u0 * current * layerRadius**2 / 2 / (((z - windingCenter)**2 + layerRadius**2)**1.5)
        B_eff = abs(B) + M * self.trapAcc * (z - z[np.argmin(abs(B))])/ (ub * mj * gj)
        return z,B,B_eff

        
    def plotField1D(self, pos, B, B_eff):
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
            ax.plot(pos*1e3, abs(B), '*', label = 'lab frame')
            ax.plot(pos*1e3, B_eff, 'r--', label = 'co-moving frame')
            ax.set_xlabel('Position[mm]')
            ax.set_ylabel('Magnetic field[T]')
            ax.set_title('Magnetic field in coil #' + str(self.trapNum))
            ax.legend()
            plt.show()
    

    
