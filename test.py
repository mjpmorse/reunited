from physics import getUsTogether
from project import noiseOn,begin

noiseOn()
begin('foo')
getUsTogether()


from numpy import loadtxt
from pylab import plot,axes,show,xlabel,ylabel,title
data = loadtxt('foo.dat')
plot(data[:,1],data[:,2],'b-')
plot(data[:,5],data[:,6],'m-')
axes().set_aspect('equal', 'datalim')
xlabel('x'),ylabel('y')
title('x,y projection of trajectories')
show()
