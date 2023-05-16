################################################################
# Computes local sideral time and transforms different coordinates
###############################################################
import math
import ephem
import datetime

###############################################################


def get_lst(lon, epoch_time):
    '''Gives local sideral time'''
    station = ephem.Observer()
    date_time = datetime.datetime.fromtimestamp(epoch_time)
    station.date = date_time
    time = station.sidereal_time()+lon
    return time

###############################################################


def get_alt(dec, rad, lat, time):
    '''Gives alt from topocentric equatorial coordinates'''
    H = time-rad
    alt = math.asin(math.sin(dec)*math.sin(lat) +
                    math.cos(lat)*math.cos(dec)*math.cos(H))
    return alt

###############################################################


def get_az(alt, dec, rad, lat, time):
    '''Gives az from topocentric equatorial coordinates'''
    H = time-rad
    az = math.acos(
        ((math.sin(dec)-math.sin(lat)*math.sin(alt))/(math.cos(lat)*math.cos(alt))))
    if (math.sin(H) > 0):
        az = 2*math.pi-az
    return az
