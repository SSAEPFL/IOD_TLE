################################################################
# Computing the orbital elements from r,v (Gibbs method)
# according to book: Orbital mechanics for engineering students
# by Howard Curtis
###############################################################
import math
import numpy as np
import constants
###############################################################


class TLE:
    def __init__(self, input):
        self.r = input[0]
        self.v = input[1]
        self.t = input[2]

        # additional TLE components
        self.SCN = self.getSCN()  # Satellite Catalog Number
        self.Class = self.getClass()  # classification of the satellite U = unclassified
        # International Designator = last two digits of launch year
        self.intDesY = self.getIntDesY()
        # International Designator = launch number of the year
        self.intDesNY = self.getIntDesNY()
        self.intDesP = self.getIntDesP()  # International Designator = piece of the launch
        self.epochY = self.getEpochY()  # last two digits of epoch year
        self.epochD = self.getEpochD()  # day and fraction of the day
        self.derMM = self.getDerMM()  # first derivative of the mean motion
        self.der2MM = self.getDer2MM()  # second derivative of the mean motion
        self.Bstar = self.getBstar()  # Drag term, radiation pressure coeff
        self.EphType = self.getEphType()  # Ephemeris type, always 0
        self.ESN = self.getESN()  # Element set number
        self.checksum = self.getChecksum()  # Checksum
        self.rev = self.getRev()  # revolution number at epoch

    def rad2deg(self, arg):
        argDeg = arg*360/(2*constants.pi)
        if argDeg < 0:
            argDeg += 360
        return argDeg

    def values(self):
        h = np.cross(self.r, self.v)
        z = [0, 0, 1]
        N = np.cross(z, h)
        rn = np.linalg.norm(self.r)
        vn = np.linalg.norm(self.v)
        vr = np.dot(self.r, self.v)/rn

        # inclination
        i = self.rad2deg(np.arccos(h[2]/np.linalg.norm(h)))

        # right ascension
        nxn = N[0]/np.linalg.norm(N)
        Omega = self.rad2deg(np.arccos(nxn))
        if N[1] < 0:
            Omega = abs(360 - Omega)

        # eccentricity - new variant
        ev = 1/constants.mu * np.cross(self.v, h) - 1/rn * self.r
        e = np.linalg.norm(ev)
        if e > 1:
            raise ValueError("Eccentricity is greater than 1")

        # argument of perigee
        omega = self.rad2deg(np.arccos(np.dot(N, ev)/(e*np.linalg.norm(N))))
        if ev[2] < 0:
            omega = 360-omega

        # true anomaly
        theta = np.arccos(np.dot(ev, self.r)/(e*rn))
        if vr < 0:
            theta = 2*constants.pi - theta

        # mean anomaly
        E = np.arcsin((np.sqrt(1-e**2)*np.sin(theta)) / (e*np.cos(theta) + 1))
        M = self.rad2deg(E - e*np.sin(E))

        # mean motion
        a = np.linalg.norm(h)**2/(constants.mu*(1-e**2))
        n = np.sqrt(constants.mu/a**3) * 3600*24/(2*constants.pi)

        return i, Omega, e, omega, M, n

    def getSCN(self):
        SCN = 0
        return f"{SCN:05d}"

    def getClass(self):
        Class = 'U'
        return Class

    def getIntDesY(self):
        IntDesY = 0
        return f"{IntDesY:02d}"

    def getIntDesNY(self):
        IntDesNY = 0
        return f"{IntDesNY:03d}"

    def getIntDesP(self):
        IntDesP = 0
        return f"{IntDesP:03d}"

    def getEpochY(self):
        time = self.t
        years = int(time / (60*60*24*365)) + 1970
        if years < 2000:
            EpochY = years - 1900
        else:
            EpochY = years - 2000
        return f"{EpochY:02d}"

    def getEpochD(self):
        time = self.t
        years = int(time / (60*60*24*365))
        EpochD = (time - (years * 60*60*24*365)) / (60*60*24)
        retprov = "{:11.8F}".format(EpochD)
        return retprov

    def getDerMM(self):
        DerMM = 0
        return "{:.7f}".format(DerMM)

    def getDer2MM(self):
        Der2MM = 0
        return f"{Der2MM:010d}"

    def getBstar(self):
        inp = 0
        #exp = fexp(inp)
        #first = int(inp * 10**(-exp+4))
        #second = exp + 1
        # return (str(first) + str(second))
        return f"{inp:08d}"

    def getEphType(self):
        EphType = '0'
        return EphType

    def getESN(self):
        ESN = 0
        return f"{ESN:04d}"

    def getChecksum(self):
        orbnum = 0
        return (orbnum % 10)

    def getRev(self):
        rev = 0
        return 0

    def list_TLE(self):
        i, Omega, e, omega, M, n = self.values()
        return [e, n, i, Omega, omega, M, self.t]

    def TLE_format(self):
        i, Omega, e, omega, M, n = self.values()

        line1 = str(1) + ' ' + self.SCN + self.Class + ' ' + self.intDesY \
            + self.intDesNY + self.intDesP + ' ' + self.epochY + self.epochD \
            + ' ' + self.derMM + ' ' + self.der2MM + ' ' + self.Bstar + ' ' \
            + self.EphType + ' ' + self.ESN + str(self.checksum)

        # inclination
        str_i = f"{i:.4f}"
        if i < 100:
            str_i = '0' + str_i
        # right ascension
        str_O = "{:8.4F}".format(Omega)
        # eccentricity
        string_e = "{:9.7F}".format(e)
        str_e = string_e[2:]
        # argument of perigee
        str_o = "{:8.4F}".format(omega)
        # mean anomaly
        str_M = "{:8.4F}".format(M)
        # mean motion
        str_n = "{:11.8F}".format(n)
        # revolution
        str_rev = "{:05d}".format(self.rev)

        line2 = str(2) + ' ' + self.SCN + ' ' + str_i + ' ' + str_O + ' '\
            + str_e + ' ' + str_o + ' ' + str_M + ' ' + str_n + str_rev \
            + str(self.checksum)

        f = open("TLE.out", "w")
        f.write(line1 + "\n" + line2)
        f.close
