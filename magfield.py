# magfield.py --- 
# 
# Filename: magfield.py
# Description: 
# 
# Author:    Yu Lu
# Email:     yulu@utexas.edu
# Github:    https://github.com/SuperYuLu 
# 
# Created: Tue Apr  2 13:49:21 2019 (-0500)
# Version: 
# Last-Updated: Wed Apr  3 21:15:27 2019 (-0500)
#           By: Yu Lu
#     Update #: 159
# 

import os
import pandas as pd
import numpy as np


class LocalField:
    """
    Magnetic field opteration class
    """
    
    def __init__(self, field, origin = (0, 0), resolution = 0.001):
        self.field = field
        self.origin = origin
        self.resolution = resolution

    @property
    def gradient(self):
        if not hasattr(self, '_gradient'):
            self._gradient = np.gradient(self.calc_absfield())
            return self._gradient
        else:
            return self._gradient
        
    def neighbor_add(self, other):
        # other: another LocalField object
        rows, cols, chann = self.field.shape # r-axis, z-axis, B-Field channels 
        newField = np.zeros([rows, cols + (cols - 1) / 2, chann])
        newField[:, :cols, :] = self.field
        newField[:, (cols - 1)/2:, :] += other.field
        return LocalField(newField, origin = self.origin, resolution = self.resolution)

    def calc_absfield(self):
        # sqrt(Bx^2 + By^2 + Bz^2)
        absfield = np.linalg.norm(self.field, 2, axis = -1)
        return absfield
                
    def get_gradient(self, xyz):
        """
        xyz: (x, y, z)
        return: (dB_r, dB_z)
        """
        x, y, z = xyz
        r = np.sqrt(x**2 + y**2)
        rIdx = np.rint((r - self.origin[0]) / self.resolution).astype(np.int32)
        zIdx = np.rint((z - self.origin[1]) / self.resolution).astype(np.int32)
        return self.gradient[0][rIdx, zIdx], self.gradient[1][rIdx, zIdx]

    
class MagField:
    """
    Magnetic field data loading and caching class 
    handles field symmetry 
    """
    def __init__(self, resolution = None, source_file = None, current = -1):
        self.source_file = source_file  # souce file for saved 2d field data
        self.resolution = resolution
        if current < 0:
            print("[!] Current cannot set to  negative !")
        self.current = current
                
    def load(self):
        if not os.path.exists(self.source_file):
            print("[!] Provided file not found !")
            return None
        
        print("[*] Loading magnetic field data ...")
        data = pd.read_csv(self.source_file)
        return data

    @property
    def cache(self):
        """
        field cache in shape [dim_r, dim_z, 3] 
        Channels are: Bx, By(=0), Bz
        """
        if not hasattr(self, '_cache'):
            data = self.load()
            
            if data is None:
                print("[!] Cannot build cache Magnetic field data not loaded !")
                return None
            
            print("[*] Caching magnetic field data ...")

            if data.z.unique().size % 2 == 0: # make sure there are odd number of z coords, otherwise drop last one
                print("[#] Found even number of field mesh along z, forcing symmetric...")
                data.drop(data.z.idxmax(), inplace = True) # double check this line 
                assert(data.z.unique().size %2 == 1)
            
            self.leftBound = data.z.min()
            self.rightBound = data.z.max()

            zs, rs = data.z.unique().size, data.x.unique().size
            cache = np.zeros([rs, zs, 3])  # Wow, this is like RGB image, channels are Bx, By, Bz

            rMin = data.x.min()
            zMin = data.z.min()
            assert(abs(rMin) < 1e-6)

            rows, cols = data.shape
            for i in range(rows):
                rIdx = int(np.rint((data.x[i] - rMin)/self.resolution))
                zIdx = int(np.rint((data.z[i] - zMin)/self.resolution))
                field = [data.Bx[i], data.By[i], data.Bz[i]] # may consider set By to 0
                cache[rIdx, zIdx, :] = field
                
            self._cache = cache
            return self._cache
        else:
            return self._cache
    
        
if __name__ == '__main__':
    mf = MagField(resolution = 5e-4, source_file = './x_y_z_Bx_By_Bz.csv', current = 400)
    field = LocalField(mf.cache, origin = (0, 0), resolution = mf.resolution)
    print(field.get_gradient((0,0,0)))
    
