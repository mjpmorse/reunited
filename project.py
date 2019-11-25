# project.py
# JR 11/15/2014

# You may import only the following functions:
#	begin, getEarthRadius, getSpeedOfLight, 
#	getGPS1Data   , walkPerson1, turnPerson1,  
#	getGPS2Data   , walkPerson2, turnPerson2,  
#	getGPSBothData, walkPersons,           
#	noiseOn, noiseOff, end
# You may not declare global any of my variables.

from numpy import *

noUncertainty = False

def noiseOff():
	global noUncertainty
	noUncertainty = True
	print "NOTE: begin() must be (re-)called after calling noiseOff()"
		
def noiseOn():
	global noUncertainty
	noUncertainty = False
	print "NOTE: begin() must be (re-)called after calling noiseOn()"
		
def begin(myruncode):
	global base
	global position1, heading1
	global position2, heading2
	global turnuncertainty1, stuncertainty, rtuncertainty, steplength1
	global turnuncertainty2,                               steplength2
	global steplengthuncertainty1, turnpermeter1
	global steplengthuncertainty2, turnpermeter2
	global basetol, ABCs, time, nGPS, logfile
	global c, earthr, Rsat
	global noUncertainty

#	noUncertainty = True # SET THIS TO True TO TURN OFF RANDOM PERTURBATIONS
	basetol = 10
	c =	 299792.4580*1000
	earthr = 6356.8171*1000
	Rsat =  20200.0000*1000
	if noUncertainty:
		print "Uncertainty turned off for debugging"
		stuncertainty = 0  # satellite clock error bound
		rtuncertainty = 0  # receiver clock error bound
		turnuncertainty1 = 0
		turnuncertainty2 = 0
		#minturnradius = 500
		turnpermeter1 = 0
		turnpermeter2 = 0
		steplengthuncertainty1 = 0
		steplengthuncertainty2 = 0
		howfaraway = 3000.
		heading1 = pi/4.
		heading2 = 0.
	else:
		stuncertainty = 1.0e-8  # satellite clock error bound
		rtuncertainty = 1.0e-4  # receiver clock error bound
		turnuncertainty1 = 2*pi/100
		turnuncertainty2 = 2*pi/100
		minturnradius = 500.
		turnpermeter1 = .005 + (2.*random.rand()-1.)*(1/minturnradius)
		turnpermeter2 = .005 + (2.*random.rand()-1.)*(1/minturnradius)
		steplengthuncertainty1 = 0.1
		steplengthuncertainty2 = 0.08
		howfaraway = 3124. 
		heading1 = random.rand()*2*pi
		heading2 = random.rand()*2*pi
	steplength1 = .95 
	steplength2 = .85 
	base = cart(earthr, 88.0123*pi/180,-123.2745*pi/180)
	position1 = base
	position2 = base
	position1 = move(position1,tangent(position1,heading1),howfaraway)
	time = 0.
	nGPS = 0
	logfile = open( myruncode+'.dat','w')
#	print >> logfile,'time latitude(degrees) longitude(degrees)  heading(degrees)\n')
	[rho,lat,lon] = sph(position1)
#	print >> logfile, '{0:12.1f} {1:24.15f} {2:24.15f} {3:24.15f}'.format(time,  position[0]-base[0], position[1]-base[1], heading*180/pi )
	print >> logfile, time,  position1[0]-base[0], position1[1]-base[1], position1[2]-base[2], heading1*180/pi , \
                             position2[0]-base[0], position2[1]-base[1], position2[2]-base[2], heading2*180/pi 

def end():
	global logfile
	logfile.close()
	

def getEarthRadius():	# returns Earth radius in meters
	return 6356.8171*1000  

def getSpeedOfLight():  # returns speed of light in meters/sec
	return 299792.4580*1000

def getBaseLatLongDegrees():
	lat_degrees  =   88.0123
	long_degrees = -123.2745
	return lat_degrees,long_degrees

def getGPSBothData(): 
	return (getGPSData(position1,timeincrement=0),getGPSData(position2))

def getGPS1Data(): 
	return getGPSData(position1)

def getGPS2Data(): 
	return getGPSData(position2)

def turnPerson1(angle):  # heading is the current angle w.r.t north
	global heading1, turnuncertainty1
	heading1 = heading1 + angle         + turnuncertainty1*(2.*random.rand()-1.)

def turnPerson2(angle):  # heading is the current angle w.r.t north
	global heading2, turnuncertainty2
	heading2 = heading2 + angle         + turnuncertainty2*(2.*random.rand()-1.)

def walkPersons(stepsPerson1,stepsPerson2):
	global position1, heading1, steplength1, logfile
	global steplengthuncertainty1, turnpermeter, base, basetol, ABCs, time, nGPS
	global position2, heading2, steplength2, logfile
	global steplengthuncertainty2, turnpermeter, base, basetol, ABCs, time, nGPS
	step1 = steplength1*(1+steplengthuncertainty1*(2.*random.rand()-1.))
	step2 = steplength2*(1+steplengthuncertainty2*(2.*random.rand()-1.))
	numsteps1 = int(stepsPerson1)
	numsteps2 = int(stepsPerson2)
	maxsteps  = max(numsteps1,numsteps2)
	for i in range(maxsteps):
		time += 1
		if i<numsteps1:
			position1 = move(position1,tangent(position1,heading1),step1)
			heading1 = heading1 + turnpermeter1*step1
		if i<numsteps2:
			position2 = move(position2,tangent(position2,heading2),step2)
			heading2 = heading2 + turnpermeter2*step2

		#[rho,lat,lon] = sph(position1)
		#print 'heading = ',heading
		print >> logfile, time,  position1[0]-base[0], position1[1]-base[1], position1[2]-base[2], heading1*180/pi , \
                               position2[0]-base[0], position2[1]-base[1], position2[2]-base[2], heading2*180/pi 
		#rh,th,ph=sph(position1) 
		#print time, rh,th*180/pi,ph*180/pi
		if(linalg.norm(position1-position2)<=basetol):
			print 'Congratulations! You have found your friend.'
			print 'It took you ', time/(60.*60.),' hours.'
			print 'You used your GPS receiver ', nGPS,' times.\n'
			quit()
	
			
def walkPerson1(steps):
	global position1, heading1, steplength1, logfile
	global steplengthuncertainty1, turnpermeter, base, basetol, ABCs, time, nGPS
	step = steplength1*(1+steplengthuncertainty1*(2.*random.rand()-1.))
	numsteps = int(steps)
	for i in range(numsteps):
		time += 1
		position1 = move(position1,tangent(position1,heading1),step)
		heading1 = heading1 + turnpermeter1*step
		#[rho,lat,lon] = sph(position1)
		#print 'heading = ',heading
		print >> logfile, time,  position1[0]-base[0], position1[1]-base[1], position1[2]-base[2], heading1*180/pi , \
                               position2[0]-base[0], position2[1]-base[1], position2[2]-base[2], heading2*180/pi 
		#rh,th,ph=sph(position1) 
		#print time, rh,th*180/pi,ph*180/pi
		if(linalg.norm(position1-position2)<=basetol):
			print 'Congratulations! You have found your friend.'
			print 'It took you ', time/(60.*60.),' hours.'
			print 'You used your GPS receiver ', nGPS,' times.\n'
			quit()
			
def walkPerson2(steps):
	global position2, heading2, steplength2, logfile
	global steplengthuncertainty2, turnpermeter, base, basetol, ABCs, time, nGPS
	step = steplength2*(1+steplengthuncertainty2*(2.*random.rand()-1.))
	numsteps = int(steps)
	for i in range(numsteps):
		time += 1
		position2 = move(position2,tangent(position2,heading2),step)
		heading2 = heading2 + turnpermeter2*step
		#[rho,lat,lon] = sph(position1)
		#print 'heading = ',heading
		print >> logfile, time,  position1[0]-base[0], position1[1]-base[1], position1[2]-base[2], heading1*180/pi , \
                               position2[0]-base[0], position2[1]-base[1], position2[2]-base[2], heading2*180/pi 
		#rh,th,ph=sph(position2) 
		#print time, rh,th*180/pi,ph*180/pi
		if(linalg.norm(position1-position2)<=basetol):
			print 'Congratulations! You have found your friend.'
			print 'It took you ', time/(60.*60.),' hours.'
			print 'You used your GPS receiver ', nGPS,' times.\n'
			quit()
			

####################################################################################
#  DO NOT CALL ANY FUNCTIONS BELOW EXCEPT POSSIBLY FOR DEBUGGING PURPOSES

def getGPSData(pos,timeincrement=5*60.): # GPS consultation takes 5 minutes
	global position, stuncertainty, rtuncertainty, time, ABCs, c, nGPS
	global nplanes, nperplane
	nGPS = nGPS + 1
	genSatelliteLocations(time)
	visible = areVisible()
	numvisible = sum(visible)
	ABCt = zeros((numvisible,4))
	d = rtuncertainty*(2*random.rand(1)-1)
	k=0
	for i in range(nperplane*nplanes):
		if visible[i]:
			ABC = ABCs[i,:]
			ABCt[k,0:3]=ABC
			approx_transit_time = linalg.norm(pos-ABC)/c + stuncertainty*(2*random.rand(1)-1) + d
			ABCt[k,3] = approx_transit_time
			k = k+1
		m,n = ABCt.shape
	print m, 'satellites visible'
	time = time + timeincrement  
#	print "For debugging, position is:",position
	return ABCt

def move(r,t,d):
 # assume r is radial vector, t is a vector orthogonal to r
  rlength = linalg.norm(r)
  tu = t/linalg.norm(t)
  alpha = d/rlength
  rnew = cos(alpha) * r + sin(alpha) * rlength * tu
  # this should still have length = Earth radius, but can drift due to round-off
  # so we renormalize here (11/18/11):
  rnew *= getEarthRadius()/linalg.norm(rnew)
  return rnew

def turnt(r,t,delta): # rotates t in plane normal to r
 # relies on t being normal to r
	ru = r/linalg.norm(r)
	tu = t/linalg.norm(t)
	n = cross(tu,ru)
	newt = cos(delta)*tu + sin(delta)*n
	return newt

def cross(u,v):	# vector cross product
	w=zeros(3)
#	print 'cross...',u,v,w
	w[0] = u[1]*v[2] - u[2]*v[1]
	w[1] = u[2]*v[0] - u[0]*v[2]
	w[2] = u[0]*v[1] - u[1]*v[0]
	return w

def sph(r):
  # phi is latitude
  # theta is longitude
	x = r[0]
	y = r[1]
	z = r[2]
	rho = linalg.norm(r)
	theta = arctan2(y,x)
	phi = arcsin(z/rho)
	return rho,theta,phi

def cart(rho,phi,theta):
	x = rho*cos(phi)*cos(theta)
	y = rho*cos(phi)*sin(theta)
	z = rho*sin(phi)
	r = array([x,y,z])
	return r
	
def northward(r):
	[rho,phi,theta]=sph(r)  
	xp = -sin(phi)*cos(theta)
	yp = -sin(phi)*sin(theta)
	zp =  cos(phi)
	t = array([xp,yp,zp])
	return t

def bearing(r,t):
	north=northward(r)
	unitt=t/linalg.norm(t)	
	beta = arccos(north*unitt)
	return beta
  
def tangent(r,bearing):
	t = northward(r)
	t = turnt(r,t,bearing)
	return t

def rotator(a):
	R = array([ [ cos(a), sin(a), 0],
					[-sin(a), cos(a), 0],
					[ 0	  , 0	  ,   1]])
	return R
	
def genSatelliteLocations(timeInSeconds):	  # time in seconds
	global ABCs
	global nplanes, nperplane, Rsat
	dayInSeconds = 60*60*24
	earthOmega = 2*pi/dayInSeconds
	satOmega = 2*earthOmega
	nplanes = 7
	nperplane = 4
	ABCs = zeros((nplanes*nperplane,3))
	for i in range(nplanes):  # set up coordinate axes on each of the orbit planes
		normal = cart(1,2*pi*i/nplanes - 2*pi*timeInSeconds/(24*60*60),1)
		# get a basis vector for orbit plane normal to earth axis
		e1 = cross([0.,0.,1.],normal)
		e1 = e1/linalg.norm(e1)
		# and the basis vector normal to that
		e2 = cross(normal,e1)
		e2 = e2/linalg.norm(e2)
		k = nperplane*(i-1)
		theta0 = i # arbitrarily set initial phase to i radians
		# nperplane satellites on each orbit, evenly spaced
		dtheta = 2*pi/nperplane
		for kk in range(nperplane):
			theta = theta0 + satOmega*timeInSeconds + kk*dtheta
			ABC = Rsat*dot( rotator(earthOmega*timeInSeconds), cos(theta)*e1+sin(theta)*e2 )
#			print 'ABCs[k+kk,:] = ',ABCs[k+kk,:]
#			print 'ABC = ', ABC
			ABCs[k+kk,:] = ABC
			k+kk
		ABCs

def areVisible():
	global ABCs
	global base
	global nplanes, nperplane
	visible = zeros(nplanes*nperplane,dtype=bool)
	for i in range(nplanes*nperplane):
		u = base
		v = ABCs[i,:]-base
		if dot(u,v) > 0:
			visible[i] = True
	return visible
	
def getPosition1():
	global position1
	return position1
def getPosition2():
	global position2
	return position2
