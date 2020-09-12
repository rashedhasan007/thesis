import numpy as np
import math

class BC:
    def fixed(nelx,nely):
        # BC's and support
        dofs=np.arange(2*(nelx+1)*(nely+1))
        fixed=np.union1d(np.array([80]),dofs[8050:8060:1])
        return fixed
