# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 00:40:46 2020

@author: Rashidul hasan (student id-1512027)
depertmant of naval architucture and marine engineering 
Bangladesh university of engineering and technology

By using this moddule we can see our desiarbale design which is created by using design module

"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt

class design_view:
    def __init__(self,x,nelx,nely):
        self.x=x
        self.nelx=nelx
        self.nely=nely
        #x=volfrac * np.ones((nely*nelx),dtype=float)
        xPhys=x.copy()
        v=-xPhys.reshape((nelx,nely)).T
        plt.ion() # Ensure that redrawing is possible
        fig,ax = plt.subplots()
        im = ax.imshow(v, cmap='gray',\
        interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
        fig.show()