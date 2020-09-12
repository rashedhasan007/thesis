from Design import design
from view import design_view
from main import optimizer
import time
from boundary import BC
from loads import load
s=design(160,80)
c=s.solid(.4)
nelx=160
nely=80
#c=s.box([100,20],60,40,.4,False,[0])
c=s.box([100,30],60,20,.4,False,[0]) 
c=s.circule(10,[200,80],.4,True,c)  
#time.sleep(5)
design_view(c,nelx,nely)
load1=load.load1(nelx,nely)
s=BC.fixed(nelx,nely)
optimizer(1e-9,1.0,.4,nelx,nely,c,load1,s)