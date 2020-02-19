"""
Convert mean anomaly to true anomaly
Date: 2020.01.20
Author: Yuri Shimane

Function converts mean anomlay to true anomaly
"""
import numpy
from scipy import optimize

def meananom2trueanom(m0,ecc):
    """function converts mean anomaly to true anomaly
    Args:
        m0 (float): mean anomaly in radians
        e (float): eccentricity
    Returns:
        (float): true anomaly in radians
    """
    # elliptical case
    if ecc>0 and ecc<1:
        def zeroEfunc(E,m0,ecc):
            fval = E - ecc*np.sin(E) - m0
            return fval
        def zeroEprime(E,ecc):
            fpval = 1 - ecc*np.cos(E)
            return fpval
        
        # initial guess of eccentric anomaly
        E0 = m0
        # Newton-raphson method to find E with initial guess E)
        E1 = optimize.newton(zeroEfunc, args=(m0,ecc), x0=E0)#, fprime=zeroEprime)
        
        theta = 2*np.arctan(np.sqrt((1+ecc)/(1-ecc))*np.tan(E1/2))
        
    else:
        theta = NaN
        
    return theta


    