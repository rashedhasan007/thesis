import numpy as np
import math
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt ;import cvxopt.cholmod

class optimizer:
    def __init__(self,Emin,Emax,volfrac,nelx,nely,x,load,B_condition):
        ndof = 2*(nelx+1)*(nely+1)

        # Allocate design variables (as array), initialize and allocate sens.
        
        xold=x.copy()
        xPhys=x.copy()
        
        g=0 # must be initialized to use the NGuyen/Paulino OC approach
        dc=np.zeros((nely,nelx), dtype=float)
        #element stiffness matrix
        #element stiffness matrix
        def lk():
        	E=1
        	nu=0.3
        	k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
        	KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        	[k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        	[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        	[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        	[k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        	[k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        	[k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        	[k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
        	return (KE)
        # FE: Build the index vectors for the for coo matrix format.
        KE=lk()
        edofMat=np.zeros((nelx*nely,8),dtype=int)
        print(edofMat)
        for elx in range(nelx):
        	for ely in range(nely):
        		el = ely+elx*nely
        		n1=(nely+1)*elx+ely
        		n2=(nely+1)*(elx+1)+ely
        		edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
        # Construct the index pointers for the coo format
        iK = np.kron(edofMat,np.ones((8,1))).flatten()
        jK = np.kron(edofMat,np.ones((1,8))).flatten()
        rmin=5.4

        # Filter: Build (and assemble) the index+data vectors for the coo matrix format
        nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
        iH = np.zeros(nfilter)
        jH = np.zeros(nfilter)
        sH = np.zeros(nfilter)
        cc=0
        for i in range(nelx):
        	for j in range(nely):
        		row=i*nely+j
        		kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
        		kk2=int(np.minimum(i+np.ceil(rmin),nelx))
        		ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
        		ll2=int(np.minimum(j+np.ceil(rmin),nely))
        		for k in range(kk1,kk2):
        			for l in range(ll1,ll2):
        				col=k*nely+l
        				fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
        				iH[cc]=row
        				jH[cc]=col
        				sH[cc]=np.maximum(0.0,fac)
        				cc=cc+1
                
                
        # Finalize assembly and convert to csc format
        H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()	
        Hs=H.sum(1)
        # BC's and support
        dofs=np.arange(2*(nelx+1)*(nely+1))
        #fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
        #fixed=np.union1d(dofs[318:2*(160+1)],np.array([2*(160+1)*(160+1)-1]))
        #fixed=dofs[0:2*(160+1)]
        fixed=B_condition
        free=np.setdiff1d(dofs,fixed)
        
        # Solution and RHS vectors
        f=np.zeros((ndof,1))
        u=np.zeros((ndof,1))
        # Set load
        f=load
        #f[320,0]=50
        #f[51841,0]=-50
        #f[25921]=-50
        #f[51601,0]=50
        #f[51681,0]=-50
        # Initialize plot and plot the initial design
        v=-xPhys.reshape((nelx,nely)).T
        plt.ion() # Ensure that redrawing is possible
        fig,ax = plt.subplots()
        im = ax.imshow(v, cmap='gray',\
        interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
        fig.show()
        print(-xPhys)
        
        def oc(nelx,nely,x,volfrac,dc,dv,g):
        	l1=0
        	l2=1e9
        	move=0.2
        	# reshape to perform vector operations
        	xnew=np.zeros(nelx*nely)
        
        	while (l2-l1)/(l1+l2)>1e-3:
        		lmid=0.5*(l2+l1)
        		xnew[:]= np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
        		gt=g+np.sum((dv*(xnew-x)))
        		if gt>0 :
        			l1=lmid
        		else:
        			l2=lmid
        	return (xnew,gt)
            
        def deleterowcol(A, delrow, delcol):
        	# Assumes that matrix is in symmetric csc form !
        	m = A.shape[0]
        	keep = np.delete (np.arange(0, m), delrow)
        	A = A[keep, :]
        	keep = np.delete (np.arange(0, m), delcol)
        	A = A[:, keep]
        	return A
        
        def sen():
            if ft==0:
                dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
            elif ft==1:
                dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
                dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
        dv = np.ones(nely*nelx)
        dc = np.ones(nely*nelx)
        ce = np.ones(nely*nelx)
        change=1
        loop=0
        penal=3.0
        ft=1
        
        def cal(g):
            # Setup and solve FE problem
        	sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
        	K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
        		# Remove constrained dofs from matrix
            # Remove constrained dofs from matrix and convert to coo
        	K = deleterowcol(K,fixed,fixed).tocoo()
        		# Solve system 
        	K = cvxopt.spmatrix(K.data,K.row.astype(np.int),K.col.astype(np.int))
        	B = cvxopt.matrix(f[free,0])
        	cvxopt.cholmod.linsolve(K,B)
        	u[free,0]=np.array(B)[:,0]
            # Objective and sensitivity
        	ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE) * u[edofMat].reshape(nelx*nely,8) ).sum(1)
        	obj=( (Emin+xPhys**penal*(Emax-Emin))*ce ).sum()
        	dc[:]=(-penal*xPhys**(penal-1)*(Emax-Emin))*ce
        
        	dv[:] = np.ones(nely*nelx)
            # Sensitivity filtering:
        	if ft==0:
        		dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
        	elif ft==1:
        		dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
        		dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
        
        		# Optimality criteria
        	xold[:]=x
        	(x[:],g)=oc(nelx,nely,x,volfrac,dc,dv,g)
        
        		# Filter design variables
        	if ft==0:   xPhys[:]=x
        	elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]
        	    
        		# Compute the change by the inf. norm
        	change=np.linalg.norm(x.reshape(nelx*nely,1)-xold.reshape(nelx*nely,1),np.inf)
        	return K,sK,obj,change
            
        loop=0
        while True:
            calculation=cal(g)
            change=calculation[3]
            loop=loop+1
            print(calculation[2],loop)
            im.set_array(-xPhys.reshape((nelx,nely)).T)
            fig.canvas.draw()
            if loop==200 or change<=.01:
                break
        
        



        