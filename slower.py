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


class Trap:
    def __init__(self, trapIdx, slowerObj):
        # Trap id starts from 0 
        self.idx = trapIdx
        self.onTime = slowerObj.trapOnTime[trapIdx]
        self.offTime = slowerObj.trapOffTime[trapIdx]
        self.pulseTime = slowerObj.trapPulseLength[trapIdx]
        self.center = slowerObj.trapCenter[trapIdx]
        self.leftCoil = slowerObj.trapLeftCoil[trapIdx]
        self.rightCoil = slowerObj.trapRightCoil[trapIdx]

    #def previous
    #def next
            
    def __str__(self):
        return '\n'.join(["%s: %.6f"%(key, value) for key,value in self.__dict__.items()])
    
class Slower:
    
    ###########
    # GLOBALS #
    ###########
    
    ## GEOMETRY OF A SINGLE TRAP COIL
    ################################
    coilSpace = 10e-3 # Measure from the center of front/back coil
        
    ## GEOMETRY OF SLOWER
    #####################
    ### Note: trap numbers: 0, 1, 2, ..., 479
    trapSpace = coilSpace/ 2. # neighboring trap center space 
    numTraps = 480             # total number of traps 
    divTrapIdx = 179           # trap number where traps divides
    
    totalLength = (numTraps - 1) * trapSpace # Total length only counting trap centers
    stage1Length = divTrapIdx  * trapSpace    
    stage2Length = totalLength - stage1Length
        
    def __init__(self, initialV, finalV, accRatio, geoOffset = 0, timeOffset = 0):
        """
        Initialize slower geometry and dynamics settings;
        All geometry parameters are calculating from trap center;
        """
        
        ## DYNAMICS OF SLOWER
        #####################
        
        ### Configurations 
        self.initialV = initialV
        self.finalV = finalV
        self.accRatio = accRatio
        
        ### 1st trap center offset and turn on time offset
        self.geoOffset = geoOffset
        self.timeOffset = timeOffset

        ### Derived configurations 
        self.stage1Acc, self.stage2Acc = self._acceleration()
        self.middleV = np.sqrt(finalV**2 - 2 * self.stage2Acc * self.stage2Length)
        self.stage1Time, self.stage2Time, self.totalTime = self._totalTime()

        ## Properties of traps
        self.trapIdx = np.arange(0, self.numTraps).astype(np.int32)
        self.trapCenter = self._trapCenter()
        self.trapLeftCoil = self._trapLeftCoil()
        self.trapRightCoil = self._trapRightCoil()
        self.trapOnTime = self._trapOnTime()
        self.trapOffTime = self._trapOffTime()
        self.trapPulseLength = self._trapPulseLength()


    def __str__(self):

        print("\nSlower Parameters:\n================================")
        print("Num Traps:  {:<25} \t Length[m]:  {:<}".format(self.numTraps, self.totalLength),\
              "\nInitial V[m/s]:  {:<25} \t Final V[m/s]:  {:<}".format(self.initialV, self.finalV),\
              "\nStage1 Acc[m/s2]:  {:<25.2f} \t Stage2 Acc[m/s2]:  {:<.2f}".format(self.stage1Acc, self.stage2Acc),\
              "\nStage1 time[ms]:  {:<25.2f} \t Stage2 time[ms]:  {:<.2f}".format(self.stage1Time*1e3, self.stage2Time*1e3),\
              "\nTotal time[ms]:  {:<25.2f} \t Middle V[m/s]:  {:<.2f}".format(self.totalTime * 1e3, self.middleV)\
              )
        return '\n'

    
    def _acceleration(self):
        """
        Calculate the acceleration value for 1st and 2nd stage
        """
        
        a1 = (self.finalV**2 - self.initialV**2) / (2 * (self.stage1Length + self.accRatio * self.stage2Length))
        a2 = a1 * self.accRatio
        return a1, a2

    def _totalTime(self):
        """
        Calculate total time for deceleration
        """
        t1 = (self.middleV - self.initialV) / self.stage1Acc
        t2 = (self.finalV - self.middleV) / self.stage2Acc
        t = t1 + t2
        return t1, t2, t
    
    
    def _trapCenter(self):
        """
        Calculate the geometric center for each trap
        """
        return self.geoOffset + self.trapSpace * self.trapIdx

    def _trapLeftCoil(self):
        return self._trapCenter() - self.trapSpace

    def _trapRightCoil(self):
        return self._trapCenter() + self.trapSpace

    def _trapVelocity(self):
        """
        Calculate the effective moving velocity of specified trap
        """
        pos = self._trapCenter()
        velocity = []
        for i in self.trapIdx:
            if i  <= self.divTrapIdx:
                v = np.sqrt(self.initialV**2 + 2 * self.stage1Acc * (pos[i] - self.geoOffset))
            else:
                v = np.sqrt(self.middleV**2 + 2 * self.stage2Acc * (pos[i] - self.geoOffset - self.stage1Length))
            velocity.append(v)
            
        return np.array(velocity)

    def _trapCenterTime(self):
        """
        Time when atoms arrive at the trap center
        """ 
        velocity = self._trapVelocity()
        times = []
        for i in self.trapIdx:
            if i <= self.divTrapIdx:
                times.append(self.timeOffset + (velocity[i]-velocity[0])/self.stage1Acc)
            else:
                times.append(self.timeOffset
                             + (velocity[self.divTrapIdx] - velocity[0]) / self.stage1Acc
                             + (velocity[i] - velocity[self.divTrapIdx]) / self.stage2Acc)
        return np.array(times)

    def _trapOnTime(self):
        return self._trapCenterTime() - self.trapSpace / self._trapVelocity()
    
    def _trapOffTime(self):
        return self._trapCenterTime() + self.trapSpace / self._trapVelocity()

    def _trapPulseLength(self):
        return self._trapOffTime() - self._trapOnTime()

    def plot(self):
        fig, ax = plt.subplots(1,2, figsize = (10,8))
        ax[0].plot(self.trapIdx, self.trapCenter, label = 'trap center')
        ax[0].set_xlabel('trap idx')
        ax[0].set_ylabel('center [m]')
        ax[1].plot(self.trapIdx, self.trapPulseLength*1e6, label = 'trap pulse length')
        ax[1].set_xlabel('trap idx')
        ax[1].set_ylabel('time [us]')
        plt.show()
        
    

if __name__ == '__main__':
    slower = Slower(480, 50, 1, geoOffset = 1, timeOffset = 1)
    trap = Trap(1, slower)
    print(slower)
    print(trap)
    slower.plot()
    
