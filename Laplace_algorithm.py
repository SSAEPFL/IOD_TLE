################################################################
# Laplace preliminary orbit determination according to book:
# Fundamentals of astrodynamics by Bate, Mueller and White
###############################################################
import math
import numpy
import numpy as np
from StationPositionVector import *
import constants
from IOD_class import *
###############################################################


def rv(iodset, n):
    '''Calcute r_n and v_n using the Laplace algorithm'''
    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "LOCATION.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    Station_Coordinates = Data.readlines()
    Data.close()
    lat1 = float(Station_Coordinates[0])*(np.pi/180)
    lon1 = float(Station_Coordinates[1])*(np.pi/180)
    h1 = float(Station_Coordinates[2])
    lat2 = float(Station_Coordinates[3])*(np.pi/180)
    lon2 = float(Station_Coordinates[4])*(np.pi/180)
    h2 = float(Station_Coordinates[5])
    lat3 = float(Station_Coordinates[6])*(np.pi/180)
    lon3 = float(Station_Coordinates[7])*(np.pi/180)
    h3 = float(Station_Coordinates[8])

    # Extracting observation time:
    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()
    time_s1 = get_lst(lon1, t1)
    time_s2 = get_lst(lon2, t2)
    time_s3 = get_lst(lon3, t3)

    # Creates the vector pointing from the center of the earth to the station
    R1 = np.array(get_R(time_s1, lat1, h1))
    R2 = np.array(get_R(time_s2, lat2, h2))
    R3 = np.array(get_R(time_s3, lat3, h3))

    # Unit vectors pointing to satellite from station
    L1 = np.array(iodset[0].get_e(lat1, time_s1))
    L2 = np.array(iodset[1].get_e(lat2, time_s2))
    L3 = np.array(iodset[2].get_e(lat3, time_s3))

    if (n == 1):
        t = t1
    elif (n == 2):
        t = t2
    else:
        t = t3

    L = (t-t2)*(t-t3)/((t1-t2)*(t1-t3)) * L1 \
        + (t-t1)*(t-t3)/((t2-t1)*(t2-t3)) * L2 \
        + (t-t1)*(t-t2)/((t3-t1)*(t3-t2)) * L3

    Ldot = (2*t-t2-t3)/((t1-t2)*(t1-t3)) * L1 \
        + (2*t-t1-t3)/((t2-t1)*(t2-t3)) * L2 \
        + (2*t-t1-t2)/((t3-t1)*(t3-t2)) * L3

    Lddot = 2/((t1-t2)*(t1-t3)) * L1 + 2/((t2-t1)*(t2-t3)) * \
        L2 + 2/((t3-t1)*(t3-t2)) * L3

    R = np.array([R1, R2, R3])
    Rdot = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])  # ToDo
    Rddot = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])  # ToDo

    # linear interpolation of Rddot
    Rd1 = (R[1]-R[0])/(t2-t1)
    Rd2 = (R[2]-R[1])/(t3-t2)
    Rdd = 2*(Rd2-Rd1)/(t3-t1)
    for i in range(3):
        Rdot[0][i] = Rd1[i]
        Rdot[1][i] = (Rd1[i]+Rd2[i])/2
        Rdot[2][i] = Rd2[i]
        Rddot[i] = Rdd

    D = 2*np.linalg.det([[L[0], Ldot[0], Lddot[0]],
                        [L[1], Ldot[1], Lddot[1]],
                        [L[2], Ldot[2], Lddot[2]]])
    D1 = np.linalg.det([[L[0], Ldot[0], Rddot[n-1][0]],
                        [L[1], Ldot[1], Rddot[n-1][1]],
                        [L[2], Ldot[2], Rddot[n-1][2]]])
    D2 = np.linalg.det([[L[0], Ldot[0], R[n-1][0]],
                        [L[1], Ldot[1], R[n-1][1]],
                        [L[2], Ldot[2], R[n-1][2]]])
    LR = np.dot(L, R[n-1])
    C1 = 4*D1**2/(D**2) - 4*D1*LR/D + np.linalg.norm(R[n-1])**2
    C2 = 8*constants.mukm*D1*D2/(D**2) - 4*constants.mukm*D2*LR/D
    C3 = 4*constants.mukm**2*D2**2/(D**2)

    # iterative process to solve the eight order equation for r2n, the norm of r2
    rn = 300 + constants.Re/1000  # first estimation: hight of ISS
    ri = 0
    while (abs(rn - ri) > 0.001):
        ri = rn
        r6 = ri**6
        r3 = ri**3
        rn = np.power(C1*r6 + C2*r3 + C3, (1./8.))

    rho = -2*D1/D - 2*constants.mukm*D2/(D*rn**3)

    r = rho*L + R[n-1]
    #calculation of the velocity
    D3 = np.linalg.det([[L[0], Rddot[n-1][0], Lddot[0]],
                        [L[1], Rddot[n-1][1], Lddot[1]],
                        [L[2], Rddot[n-1][2], Lddot[2]]])
    D4 = np.linalg.det([[L[0], R[n-1][0], Lddot[0]],
                        [L[1], R[n-1][1], Lddot[1]],
                        [L[2], R[n-1][2], Lddot[2]]])
    rhodot = -D3/D - constants.mukm*D4/(D*(np.linalg.norm(r)**3))
    v = rhodot*L + rho*Ldot + Rdot[n-1]
    return r, v

###############################################################
