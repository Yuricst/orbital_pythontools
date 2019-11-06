"""
Lambert algorithm
Date: 2019.10.12
Author: Yuri Shimane

Script contains algorithm to solve Lambert's Problem for orbital transfers as well as 
definitions of functions used in the derivation of the solution. 
Formulation follows the derivation from chapter 5.3 in Curtis (2014) "Orbital Mechanics for Engineering Students 3/e"

"""

import numpy as np
from numpy import linalg as LA

def Stumpff_S(z):
	"""
	Stumpff function S(z)
	Args:
		z (float): universal anomaly^2/semi-major axis of transfer trajectory
	Returns:
		(float): value of Stumpff functio S(z) evaluated for input z
	"""
	if z > 0:
		S = (np.sqrt(z) - np.sin(np.sqrt(z)))/np.power(z,1.5)
	elif z == 0:
		S = 1/6
	elif z < 0:
		S = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.power(-z,1.5)

	return S


def Stumpff_C(z):
	"""
	Stumpff function C(z)
	Args:
		z (float): universal anomaly^2/semi-major axis of transfer trajectory
	Returns:
		(float): value of Stumpff functio S(z) evaluated for input z
	"""
	if z > 0:
		C = (1 - np.cos(np.sqrt(z)))/z
	elif z == 0:
		C = 1/2
	elif z < 0:
		C = (np.cosh(np.sqrt(-z)) - 1)/(-z)

	return C


def y_538(r1,r2,A,z):
	"""
	Intermediate function in Lambert problem derivation (eq.5.38 in Curtis)
	Args:
		r1 (1x3 numpy array): position vector of departure
		r2 (1x3 numpy array): position vector of arrival
		A (float): intermediate value related only to input parameters
		z (float): universal anomaly^2/semi-major axis of transfer trajectory
	Returns:
		(float): value of function evaluated for input z
	"""

	y = LA.norm(r1) + LA.norm(r2) + A*(z*Stumpff_S(z) - 1)/np.sqrt(Stumpff_C(z))

	return y


def lambert(r1,r2,tof,mu,grade=None):
	"""
	Function takes in classic parameters to Lambert problem to determine orbitalelements
	Args:
		r1 (1x3 numpy array): initial position vector of departure [km]
		r2 (1x3 numpy array): final position vector of arrival [km]
		tof (float): time of flight [s]
		mu (float): gravitational parameter [km^3/s^2]
		grade (str): trajectory orientation ('pro' for prograde or 'retro' for retrograde). If not provided, assume prograde orbit.
	Returns:
		(tuple): velocity vector at position 1 and 2 of solution trajectory to Lambert problem
	"""

	# compute dtheta [radians]
	tmp = np.cross(r1,r2)
	if grade=='retro':
		if tmp[2] < 0:
			dtheta = np.arccos(np.dot(r1,r2)/(LA.norm(r1)*LA.norm(r2)))
		else: 
			dtheta = 2*np.pi - np.arccos(np.dot(r1,r2)/(LA.norm(r1)*LA.norm(r2)))
	else:
		if tmp[2] < 0:
			dtheta = 2*np.pi - np.arccos(np.dot(r1,r2)/(LA.norm(r1)*LA.norm(r2)))
		else:
			dtheta = np.arccos(np.dot(r1,r2)/(LA.norm(r1)*LA.norm(r2)))
	
	print('dtheta: {}'.format(dtheta*360/(2*np.pi)))

	# compute input parameter A where A = sin(dtheta) * sqrt[r1*r2 / (1 - cos(dtheta))]
	A = np.sin(dtheta) * np.sqrt(LA.norm(r1)*LA.norm(r2)/(1 - np.cos(dtheta)))

	print('Value of A: {}'.format(A))

	# Newton's method to find z
	# initial guess for z
	z = 1.5                #FIXME - with guess of orbital element?
	eps = 0.00001 # define error value for Newton's method
	diff = 1
	count = 1

	while diff > eps:
		print('current iteration: {}'.format(count))
		# define functions to solve by Newton's method
		F = np.power(y_538(r1,r2,A,z)/Stumpff_C(z), 3/2) * Stumpff_S(z) + A*np.sqrt(y_538(r1,r2,A,z)) - np.sqrt(mu)*tof

		if z == 0:
			Fdot = np.sqrt(2) * np.power(y_538(r1,r2,A,0),1.5)/40 + (A/8)*(np.sqrt(y_538(r1,r2,A,0)) + A*np.sqrt(1/(2*y_538(r1,r2,A,0))))
		else:
			Fdot = np.power(y_538(r1,r2,A,z)/Stumpff_C(z), 1.5) * (((1/(2*z)) * (Stumpff_C(z) - 3*Stumpff_S(z)/(2*Stumpff_C(z)))) + 3*np.power(Stumpff_S(z),2)/(4*Stumpff_C(z))) + (A/8)*(3*Stumpff_S(z)*np.sqrt(y_538(r1,r2,A,z))/Stumpff_C(z) + A*np.sqrt(Stumpff_C(z)/y_538(r1,r2,A,z)))

		# iterate for next term
		z1 = z - F/Fdot
		print('current z: {}'.format(z1))

		# update difference term
		diff = np.abs(z - z1)

		# update current guess
		z = z1

		# update current iteration number
		count += 1

		#FIXME: include escape clause
		if count > 1000:
			print('Newton\'s method exceeded {} steps; consider changing initial guess of z for better result'.format(count))
			break


	# display orbit type
	if z > 0:
		print('Transfer trajectory is an ellipse; z = {}'.format(z))
	elif z == 0:
		print('Transfer trajectory is a parabola; z = {}'.format(z))
	elif z < 0:
		print('Transfer trajectory is a hyperbolla; z = {}'.format(z))

	# calculate Lagrange functions
	f = 1 - y_538(r1,r2,A,z)/LA.norm(r1)
	g = A*np.sqrt(y_538(r1,r2,A,z)/mu)
	gdot = 1 - y_538(r1,r2,A,z)/LA.norm(r2)

	print(f,g,gdot)

	# calculate initial and final velocity vectors
	v1 = (1/g)*(r2 - f*r1)
	v2 = (1/g)*(gdot*r2 - r1)
	print('Velocity at r1: {}, velocity at r2: {}'.format(v1,v2))

	return v1, v2




