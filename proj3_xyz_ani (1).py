import matplotlib as mpl
mpl.use("AGG")
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

#unpack and load data from outTimeSpreadaheet, partNum is for the particle numbers, stepNum is for time steps, etc.
partNum,stepNum,t,x,y,z,vt,vx,vy,vz,at,ax,ay,az=np.loadtxt("outTimeSpreadsheet.csv",delimiter='\t', unpack=True)

#total steps and total particles
MSTEP=int(max(stepNum)+1)
MPART=int(max(partNum)+1)


#create figure object instance and Axes3D object instace
fig=plt.figure()
ax=fig.gca(projection='3d')

#axes labels
plt.xlabel('x pos (?)')
plt.ylabel('y pos (?)')
ax.set_zlabel("z pos (?)")

#axis limits
lim=3
ax.axis([-lim, lim, -lim, lim])
ax.set_zlim3d([-lim,lim])

#scale correction, atrempt to remove squishing on z-axis
ax.auto_scale_xyz
ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 1, 1.37, 1])) 

#line object initialization
line, = ax.plot([0],[0],[0],'.', markersize=1)
line1, = ax.plot([0],[0],[0],'o', markersize=5)
line2, = ax.plot([0],[0],[0],'o', markersize=5)

#special particle numbers for highlighting, 0 is invader, 1 is nucleus
p1=1
p2=0

for n in range(0, MSTEP):
	#first and last elements to read from during this timestep
	first=n*MPART
	last=(first+MPART)
	
	#defining the arrays holding values to plot
	xn=list(x[first:last])
	yn=list(y[first:last])
	zn=list(z[first:last])
	
	#below: editing data of line object is computationally effecient over using plot(), plot() redraws everything, 
	#and apparently axes take a long time to deaw for some reaaon (maybe because calculating 3d to 2d projection is ecpensive?)
	
	#set data of default line object 
	line.set_data(xn,yn)
	line.set_3d_properties(zn)
	#set line data for p1 and p2
	line1.set_data(xn[p1],yn[p1])
	line1.set_3d_properties(zn[p1])
	line2.set_data(xn[p2],yn[p2])
	line2.set_3d_properties(zn[p2])
	
	#save images in images folder (make sure you have images folder first)
	fig.canvas.flush_events()
	fig.savefig("images/%04d.png" % n)
	

