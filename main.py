import Gauss_algorithm
import Laplace_algorithm
from IOD_class import *
from TLE_class import *
import math
import numpy as np
###############################################################


def main():
    # Creates three IOD
    IODset = np.array(Create_IODs())
###############################################################
# Computes speed and position of satellite with Gauss algorithm
###############################################################
    # precision
    prec = 0.0001
    # maximal iterations
    itmax = 1000
    r2, v2 = Gauss_algorithm.find_r(IODset, itmax, prec, True)
###############################################################
# Computes speed and position of satellite with Laplace algorithm
###############################################################
    r2Lapl, v2Lapl = Laplace_algorithm(IODset, 2)

###############################################################
# Computes the TLE
###############################################################
    # Corresponding observatioPositionn time
    t2 = IODset[1].get_time()
    TLEinput = [r2, v2, t2]
    # Create TLE
    tle = TLE(TLEinput)
    # Create TLE file
    tle.TLE_format()


if __name__ == "__main__":
    main()
