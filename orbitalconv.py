"""
State vector 2 Orbital Elements
Date: 2019.10.12
Author: Yuri Shimane

Script contains function to convert state vector to orbital elements

"""

import numpy as np
from numpy import linalg as LA


class sv2el_return:
    """class holds output to sv2el
    """
    def __init__(self,h,i,RAAN,e,omega,theta):
        self.h = h           # specific angular momentum
        self.i = i           # inclination
        self.RAAN = RAAN     # right ascension of ascending node
        self.e = e           # eccentricity
        self.omega = omega   # argument of periapsis
        self.theta = theta   # true anomaly

def sv2el(r,v,mu):
	"""
	Funciton calculates orbital elements from state vector.
	Args:
		r (1x3 numpy array): position vector of object [km/s]
		v (1x3 numpy array): velocity vector of object [km/s]
		mu (float): gravitational parameter [km^3/s^2]
	Returns:
		(class): class holding orbital elements
			h: specific angular momentum [km^2/s]
			i: inclination [rad]
			RAAN: right ascension of ascending node [rad]
			e: eccentricity
			omega: argument of periapsis [rad]
			theta: true anomaly [rad]
	"""
	
	# specific angular momentum vector
	h_vect = np.cross(r,v)
	h = LA.norm(h_vect) # scalar-form of specific angular momentum vector

	# inclination angle
	i = np.arccos(h_vect[2]/LA.norm(h_vect))

	# eccentricity vector
	e_vect = np.cross(v,h_vect)/mu - r/LA.norm(r)
	e = LA.norm(e_vect) # scalar-form of eccentricity

	# right ascension
	K = np.array([0,0,1])
	N = np.cross(K,h_vect)
	if N[1] > 0:  # if N_y > 0
		RAAN = np.arccos(N[0]/LA.norm(N))
	else:         # if N_y <= 0
		RAAN = 2*np.pi - np.arccos(N[0]/LA.norm(N))

	# argument of periapsis
	if e_vect[2] > 0:  # if e_z > 0 then omega < 180
		omega = np.arccos(np.dot(e_vect,N)/(LA.norm(e_vect)*LA.norm(N)))
	else:              # if e_z < 0 then omega > 180
		omega = 2*np.pi - np.arccos(np.dot(e_vect,N)/(LA.norm(e_vect)*LA.norm(N)))

	# true anomaly
	v_radial = np.dot(v,r)/LA.norm(r) # radial velocity
	if v_radial > 0:  # if v_radial > 0, theta<180 (moving away)
		theta = np.arccos(np.dot(e_vect,r)/(LA.norm(e_vect)*LA.norm(r)))
	else:             # if v_radial < 0, theta<180 (moving away)
		theta = 2*np.pi - np.arccos(np.dot(e_vect,r)/(LA.norm(e_vect)*LA.norm(r)))


	# return orbital elements as dictionary
	return sv2el_return(h,i,RAAN,e,omega,theta)
