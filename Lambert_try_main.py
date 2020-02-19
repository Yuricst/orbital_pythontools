"""
main function to run lambert_alg.py
Date: 2019.10.12
Author: Yuri Shimane
"""

import numpy as np
import lambert as lb
import orbitalconv as oc

# initial and final position vectors [km]
r_in = np.array([5000, 10000, 2100])
r_fn = np.array([-14600, 2500, 7000])

# time of flight [s]
dt = 3600

# gravitational parameter [km^3/s^2]
mu_E = 398600

# solve Lambert problem
v1, v2 = lb.lambert(r1=r_in, r2=r_fn, tof=dt, mu=mu_E, grade='pro')

# obtain orbital elements
elements_in = oc.sv2el(r_in,v1,mu_E)

# print example
print(elements_in.h)


