import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl

########## Antenna in star shape pattern for running simulations
########## Taken from Anne script - 12/02/2017 OMH

def computeAntennaPos(zenz,azz,decay,dist_toxmax,dist_fromxmax):
  # Here (x,y,z) = (Northing,Westing, Up). Alt ref at sea level.
  # Use Zhaires (theta, phi) convention: defined wrt direction of origin


  DISPLAY = 0
  #pl.ion()
  deg2rad = np.pi/180

  #fname="coord_{0}_{1}_{2}_{3}.txt".format(int(zen*10),int(az),int(decay[2]),int(dist_fromxmax))
  fname = "AntennaPositions.inp"
  file= open(fname,'w') #'a')# 

  # Angles in GRAND conventions
  zen = 180-zenz
  az = 180+azz

  az_rad = az*deg2rad
  zen_rad= zen*deg2rad

  # Direction where Earth mag field points to @ Ulastai
  #az_B = 2.66*deg2rad
  az_B = 0*deg2rad  # North = Magnetic North
  zen_B = 153.18*deg2rad  #Direction where we point to
  B = np.array([np.cos(az_B)*np.sin(zen_B),np.sin(az_B)*np.sin(zen_B),np.cos(zen_B)]) #in LOFAR coordinates

  v = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])
  vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
  vxB = vxB/np.linalg.norm(vxB)
  vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
  vxvxB = vxvxB/np.linalg.norm(vxvxB)
  print "vxB:",vxB
  print "vxvxB:",vxvxB

  plane0 = (dist_toxmax+dist_fromxmax)*np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])+decay
  print 'Center of antenna plane:',plane0
  max_ang = 3*deg2rad  # Most distant antenans are 2degs from axis
  ### to keep antenna density high eneough inside and on the cherenkov cone, 3deg opening fro xmax
  d = dist_fromxmax*np.tan(max_ang) # old: cone opens at Xmax
  step = d/20.   # 20 antennas per arm (160 in total)
  #step = 25   #AZ setup
  print 'Inside/On/close by: Antenna step = ',step, ' width of cone at ', dist_fromxmax, 'm distance to Xmax =', d , ' size of Ch.cone ',  dist_fromxmax*np.tan(1.4*deg2rad), ' former size of cone ', dist_fromxmax*np.tan(max_ang)
  
  ### add a few antenna outside that ring, to be consistent with cone selection, 3 antennas added
  d2= (dist_toxmax+dist_fromxmax)*np.tan(max_ang) # new: cone opens at decay point, like in cone selection
  step2= (d2-d)/3.
  

  if DISPLAY:
  	fig = pl.figure(1)
  	#ax = fig.add_subplot(221, projection='3d')
	ax = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)


  #for i in np.arange(1,16):
  for i in np.arange(1,21):   #AZ setup
    for j in np.arange(8):
  	xyz0 = i*step*(np.cos(j/4.0*np.pi)*vxB+np.sin(j/4.0*np.pi)*vxvxB) # distances between the antennas on one arm =25m
  	xyz = plane0+xyz0  # Translation to GRAND referential
        if DISPLAY:
	  #print "Adding antenna at ",xyz
          #ax.scatter(xyz[0],xyz[1],xyz[2])
  	  xplane = np.dot(xyz0,vxB)
  	  yplane = np.dot(xyz0,vxvxB)
  	  
	  #ax.scatter(xplane,yplane)
  	  #ax2.scatter(xyz[1],xyz[2])
	  #ax3.scatter(xyz[0],xyz[2])
	  #ax4.scatter(xyz[0],xyz[1])


  	file.write('AddAntenna {0} {1} {2}\n'.format(xyz[0],xyz[1],xyz[2]))
  	
    for i in np.arange(1,4):
      for j in np.arange(8):

        xyz0 = (d+ i*step2) *(np.cos(j/4.0*np.pi)*vxB+np.sin(j/4.0*np.pi)*vxvxB) # distances between the antennas on one arm =25m
        xyz = plane0+xyz0  # Translation to GRAND referential
        file.write('AddAntenna {0} {1} {2}\n'.format(xyz[0],xyz[1],xyz[2]))

  file.close()
  if DISPLAY:
    #ax.set_ylabel('#ax.set_xlabel('Northing [m]')
    #ax.set_ylabel('Westing [m]')
    #ax.set_zlabel('Up [m]')
    ax.set_xlabel('vxB [m]')
    ax.set_ylabel('vx(vxB) [m]')
    ax.grid(True)
    ax2.set_xlabel('y [m]')
    ax2.set_ylabel('z [m]')
    ax2.grid(True)
    ax3.set_xlabel('x [m]')
    ax3.set_ylabel('z [m]')
    ax3.grid(True)
    ax4.set_xlabel('x [m]')
    ax4.set_ylabel('y [m]')
    ax4.grid(True)

    #pl.show()
    pl.savefig("star.png")
if __name__ == '__main__':
  if np.size(sys.argv) != 8:
    print "Arguments = zen az [decay] dist_toxmax dist_fromxmax. zen = 90 deg at horizon, az = 0 deg at North (>0 westward), [decay] = injection point coords (x=SN,y=EW,z=up asl), and dist_toxmax in m from decay, dist_fromxmax in m from xmax to plane."
  else:
    theta = float(sys.argv[1])
    phi = float(sys.argv[2])
    decay = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
    dist_toxmax = float(sys.argv[6])
    dist_fromxmax = float(sys.argv[7])
    print "****Shower direction = (",theta,phi,")deg, decay point = ",decay,", distance from decay=",dist_toxmax,'+',dist_fromxmax,"m"
    computeAntennaPos(theta,phi,decay,dist_toxmax,dist_fromxmax)
