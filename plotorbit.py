"""
Orbit plot
Date: 2019.10.13
Author: Yuri Shimane

Function plots orbit in 3D from Keplerian orbital elements

"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def RotMat1(phi):
	"""
	Rotational matrix about first axis
	Args:
		phi (float): angle in radians to rotate about first axis
	Returns:
		(3x3 numpy matrix): rotational matrix about first axis
	"""

	R1 = np.matrix([[1, 0, 0], [0, np.cos(phi), np.sin(phi)], [0, -np.sin(phi), np.cos(phi)]])
	return R1

def RotMat2(phi):
	"""
	Rotational matrix about second axis
	Args:
		phi (float): angle in radians to rotate about first axis
	Returns:
		(3x3 numpy matrix): rotational matrix about first axis
	"""

	R2 = np.matrix([[np.cos(phi), 0, -np.sin(phi)], [0, 1, 0], [np.sin(phi), 0, np.cos(phi)]])
	return R2

def RotMat3(phi):
	"""
	Rotational matrix about third axis
	Args:
		phi (float): angle in radians to rotate about first axis
	Returns:
		(3x3 numpy matrix): rotational matrix about first axis
	"""

	R3 = np.matrix([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
	return R3


def plot3d(mu,h,e,i,RAAN,omega,theta_range=None,interp=100,center=None):
	"""
	Function plots orbit in 3D space given orbital parameters
	Args:
		mu (float): gravitational parameter of orbiting body [kg^3/s^2]
		h (float): semi-major axis [km^2/s]
		e (float): eccentricity vector
		i (float): inclination [rad]
		RAAN (float): right ascension of ascending node [rad]
		omega (float): argument of periapsis [rad]
		theta (1x2 array): range of true anomaly to plot orbit, given in array [rad]
		interp (int): integer number of points to interpolate the trajectory
		center (float): optional radius of center body
	Returns:
		(plotly): plotly 3D plot of trajectory
		
	"""

	if theta_range:
		# plot orbit between the ranges
		theta = np.linspace(theta_range[0],theta_range[1],num=interp)
	else:
		if e < 1:
			# plot orbit over one revolution
			theta = np.linspace(0, 2*np.pi, num=interp)
		else:
			sys.exit('For escape orbit, specify range of theta to plot the orbit')

	# initialize
	r_GECtmp = np.zeros((3,1))
	r_PF = np.zeros((2,interp))
	r_GEC = np.zeros((3,interp))

	for i in range (interp):
		# compute position vector in perifocal frame
		A = (h**2)/(mu*(1 + e*np.cos(theta[i])))
		# x-coordinate in perifocal frame
		r_GECtmp[0] = A*np.cos(theta[i])
		# y-coordinate in perifocal frame
		r_GECtmp[1] = A*np.sin(theta[i])

		# extract 2D vector
		r_PF[0,i] = r_GECtmp[0]
		r_PF[1,i] = r_GECtmp[1]

		# convert vector to GEC frame
		tmp1 = np.dot(RotMat3(-omega), r_GECtmp)   #FIXME - ISSUE HERE?
		tmp2 = np.dot(RotMat1(-i), tmp1)
		tmp3 = np.dot(RotMat3(-RAAN), tmp2)

		# store vector to array
		r_GEC[0,i] = tmp3[0]
		r_GEC[1,i] = tmp3[1]
		r_GEC[2,i] = tmp3[2]


	# plot orbit in 2D
	plt.plot(r_PF[0,:], r_PF[1,:], label='Lambert transfer', color='b')
	plt.xlabel('x [km]')
	plt.ylabel('y [km]')
	plt.title('2D orbit plot in Perifocal Frame')
	plt.axis('equal')
	plt.grid(color='k', linestyle='--', linewidth=0.1)
	plt.legend()
	plt.show()

	# plot trajectory in 3D

	fig = plt.figure(figsize=(9, 9))
	ax = fig.add_subplot(111, projection = '3d')
	ax.plot(r_GEC[0,:], r_GEC[1,:], r_GEC[2,:])
	plt.title('Orbit plot in 3D')
	plt.show()


	# return position vectors in Perifocal frame and GEC
	return r_PF, r_GEC





		








