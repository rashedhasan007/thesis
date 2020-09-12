# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 00:40:46 2020

@author: Rashidul hasan (student id-1512027)
depertmant of naval architucture and marine engineering 
Bangladesh university of engineering and technology

By using this moddule we can create our desiarbale design and it return design variable x

"""
import math
from scipy.sparse import coo_matrix
import numpy as np
class design:
    def __init__(self,nelx,nely):
        self.nelx=nelx
        self.nely=nely
    # for solid object
    def solid(self,volfrac):
        nelx=self.nelx
        nely=self.nely
        x=volfrac * np.ones((nely*nelx),dtype=float)
        return x
    #for trimming  circle in a solid object 
    def circule(self,redius,center,volfrac,init,c):
        nelx=self.nelx
        nely=self.nely
        r1=center[0]
        r2=center[1]
        if init==False:
            x=volfrac * np.ones((nely*nelx),dtype=float)
            x=np.reshape(x,(-1,nely))
        else:
            x=c
            x=np.reshape(x,(-1,nely))
        for i1 in range(nelx):
            for i2 in range(nely):
                if math.sqrt((i1-(r1/2))**2+((i2-r2/2)**2))<redius:
                    x[i1,i2]=0.001
        x=np.reshape(x,(-1,1))
        x=np.array(x).flatten()
        return x
    # for trinning box in a solid object
    def box(self,point,length,breadth,volfrac,init,c):
        nelx=self.nelx
        nely=self.nely
        if init==False:
            x=volfrac * np.ones((nely*nelx),dtype=float)
            x=np.reshape(x,(-1,nely))
        else:
            x=c
            x=np.reshape(x,(-1,nely))
        for i3 in range(point[0],point[0]+length):
            for i4 in range(point[1],point[1]+breadth):
                x[i3,i4]=.001
        x=np.reshape(x,(-1,1))
        x=np.array(x).flatten()
        return x
