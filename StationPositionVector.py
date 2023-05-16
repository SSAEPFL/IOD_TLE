################################################################
# Computing the vector pointing from earth center to the
# telescope position needed for Gauss_algorithm.py according to book:
# Orbital mechanics for engineering students by Howard Curtis
###############################################################
import datetime
import math
import numpy as np
import ephem
import constants
###############################################################


def get_R(time, lat, H):
    # read latitudes GPS coordinates
    # in km!
    R = (((constants.Re/(math.sqrt(1-(2*constants.f-constants.f*constants.f)*math.sin(lat)*math.sin(lat)))+H)*math.cos(lat)*math.cos(time))/1000,
         ((constants.Re/(math.sqrt(1-((2*constants.f-constants.f*constants.f)
          * math.sin(lat)*math.sin(lat))))+H)*math.cos(lat)*math.sin(time))/1000,
         ((constants.Re*(1-constants.f)*(1-constants.f)/(math.sqrt(1-(2*constants.f-constants.f*constants.f)*math.sin(lat)*math.sin(lat)))+H)*math.sin(lat))/1000)
    return R
