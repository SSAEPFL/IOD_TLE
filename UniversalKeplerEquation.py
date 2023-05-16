################################################################
# Solving the universal Kepler equation according to book:
# Orbital mechanics for engineering students by Howard Curtis
# needed for Gauss_algorithm.py
################################################################
import math
import numpy as np
from constants import mukm


def S(z):
    if z > 0:
        Sz = (math.sqrt(z)-math.sin(math.sqrt(z)))/math.sqrt(z)**3
        return Sz
    if z < 0:
        Sz = (math.sinh(math.sqrt(-z))-math.sqrt(-z))/math.sqrt(-z)**3
        return Sz
    if z == 0:
        return 1.0/6

###############################################################


def C(z):
    if z > 0:
        Cz = (1-math.cos(math.sqrt(z)))/z
        return Cz
    if z < 0:
        Cz = (math.cosh(math.sqrt(-z))-1)/(-z)
        return Cz
    if z == 0:
        return 1.0/2

###############################################################


def get_chi(t, r, v):
    '''Solves universal Kepler's eq for
    universal anomaly, r,v_r must be in km and km/s
    otherwise overflow error'''

    # computes norm
    r_mag = math.sqrt(np.dot(r, r))
    v_mag = math.sqrt(np.dot(v, v))
    # Reciprocal of semi major axis
    alpha = 2/r_mag-v_mag**2/mukm

    # radial component of speed
    v_rad = np.dot(v, r)/r_mag
    # Initial guess of universal anomaly
    Ai = math.sqrt(mukm)*abs(alpha)*t
    # This is the tolerance of the algorithm
    epsilon = 10**(-8)
    ratio = 2*epsilon
    while abs(ratio) > epsilon:
        zi = alpha*Ai**2
        f = r_mag*v_rad/math.sqrt(mukm)*Ai**2*C(zi)+(1-alpha*r_mag) * \
            Ai**3*S(zi)+r_mag*Ai-math.sqrt(mukm)*t
        f_prime = r_mag*v_rad/math.sqrt(mukm)*Ai*(1-alpha*Ai**2*S(zi)) + \
            (1-alpha*r_mag)*Ai**2*C(zi)+r_mag
        ratio = f/f_prime
        Ai = Ai-ratio
    return Ai
