import numpy as np
import math

class load:
    def load1(nelx,nely):
    # dofs:
        ndof = 2*(nelx+1)*(nely+1)
    
        # Solution and RHS vectors
        f=np.zeros((ndof,1))
        u=np.zeros((ndof,1))
        # Set load
        #f[51521,0]=50
        #f[51841,0]=-50
        #f[ndof-121,0]=50 ##clip
        #f[ndof-41,0]=-50 ##clip
        f[ndof-101,0]=50
        f[ndof-61,0]=-50
        print(f.shape)
        return f
