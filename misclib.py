# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

MISCELLANOUS LIBRARY:
Library of Python functions for all sort of things, mostly: conversions, georeferencing specificities, naming conventions...

LIST:

- along_track_Terra: approximates along track angle of Terra orbit at a given latitude

- latlon_to_UTM: extract UTM zone and EPSG projection code from lat/lon ; for zone edges, priority goes toward SRTM naming convention (referenced to South West corner of tile)

- extended_factor_latitude_L1A: derive a factor accounting for distorsion between EPSG:4326 and ground distance
here based on : ASTER L1A granule size of 60x60km, approximated along track angle of Terra, latitude

- datetime_to_yearfraction: convert datetime to decimal years

- SRTMGL1_naming_to_latlon: extract lat/lon from SRTMGL1 naming convention (South West corner of tile) for 1x1° degree tiles

- progress_bar: progress bar for stdout readability
"""
import sys
# import utm
import numpy as np
import math as m
from datetime import datetime as dt
import time
import pandas as pd

def latlon_to_SRTMGL1_naming(lat,lon):

    if lat<0:
        str_lat = 'S'
    else:
        str_lat = 'N'

    if lon<0:
        str_lon = 'W'
    else:
        str_lon = 'E'

    tile_name = str_lat+str(int(abs(np.floor(lat)))).zfill(2)+str_lon+str(int(abs(np.floor(lon)))).zfill(3)

    return tile_name

def SRTMGL1_naming_to_latlon(tile_name):
    if tile_name[0] == 'S' or tile_name[0] == 's':
        lat = -int(tile_name[1:3])
    elif tile_name[0] == 'N' or tile_name[0] == 'n':
        lat = int(tile_name[1:3])
    else:
        sys.exit('Could not read latitude according to SRTMGL1 naming convention.')

    if tile_name[3] == 'W' or tile_name[3] == 'w':
        lon = -int(tile_name[4:7])
    elif tile_name[3] == 'E' or tile_name[3] == 'e':
        lon = int(tile_name[4:7])
    else:
        sys.exit('Could not read longitude according to SRTMGL1 naming convention.')

    return lat, lon


def datetime_to_yearfraction(date):

    #ref:https://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    def sinceEpoch(date): #returns seconds since epoch

        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def extended_factor_latitude_L1A(latitude, along_track_angle):
    def degree_lon_to_meters(lat):  #function of latitude only

        #WGS84 ellipsoid: https://gis.stackexchange.com/questions/75528/understanding-terms-in-length-of-degree-formula/75535#75535

        a = 6378137  #principal radius in meters
        f = 298.257223563  #inverse flattening

        e2 = (2 - 1 / f) / f  #squared eccentricity

        N = a / (1 - e2 * m.sin(m.pi * lat / 180) ** 2) ** 0.5  #radius of curvature of parallel at latitude lat

        R = m.cos(m.pi * lat / 180) * N  #radius of the parallel

        dist = m.pi * R / 180  #in meters

        return dist

    def degree_lat_to_meters(lat):  #function of latitude only

        #WGS84 ellipsoid: https://gis.stackexchange.com/questions/75528/understanding-terms-in-length-of-degree-formula/75535#75535

        a = 6378137  #principal radius in meters
        f = 298.257223563  #inverse flattening

        e2 = (2 - 1 / f) / f  #squared eccentricity

        M = a * (1 - e2) / (
                1 - e2 * m.sin(m.pi * lat / 180) ** 2) ** 1.5  #radius of curvature of meridian at latitude lat

        dist = m.pi * M / 180  #in meters

        return dist

    center_scene_to_ref = 0.5  #in degrees, 0.5 if using 1x1 degree ref tiles

    extent_scene = 150 * 1000  #in m, usually 60km for ASTER but using a 30% margin
    max_extent_scene = 0.5 * extent_scene * (
                m.cos(m.pi * along_track_angle / 180) + m.sin(m.pi * along_track_angle / 180))  # in m

    ext_lon = center_scene_to_ref + max_extent_scene / degree_lon_to_meters(latitude)
    ext_lat = center_scene_to_ref + max_extent_scene / degree_lat_to_meters(latitude)

    return ext_lon, ext_lat  #extended factor in degrees necessary for each direction (to implement North & South and East & West)


def latlon_to_UTM(lat,lon):

    #utm module excludes regions south of 80°S and north of 84°N, unpractical for global vector manipulation
    # utm_all = utm.from_latlon(lat,lon)
    # utm_nb=utm_all[2]

    #utm zone from longitude without exclusions
    if -180<=lon<180:
        utm_nb=int(np.floor((lon+180)/6))+1 #lon=-180 refers to UTM zone 1 towards East (West corner convention)
    else:
        sys.exit('Longitude value is out of range.')

    if 0<=lat<90: #lat=0 refers to North (South corner convention)
        epsg='326'+str(utm_nb).zfill(2)
        utm_zone=str(utm_nb).zfill(2)+'N'
    elif -90<=lat<0:
        epsg='327'+str(utm_nb).zfill(2)
        utm_zone=str(utm_nb).zfill(2)+'S'
    else:
        sys.exit('Latitude value is out of range.')

    return epsg, utm_zone


def progress_bar(value, endvalue, bar_length=20):
    #ref: https://stackoverflow.com/questions/6169217/replace-console-output-in-python

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()


def along_track_Terra(lat):
    #ref: http://fgg-web.fgg.uni-lj.si/~/mkuhar/Pouk/SG/Seminar/Vrste_tirnic_um_Zemljinih_sat/Orbit_and_Ground_Track_of_a_Satellite-Capderou2005.pdf

    nu=14.57 #nodal period

    phi=lat*m.pi/180 #latitude

    i=98.2098*m.pi/180 #inclination angle

    num=m.cos(i)-(1/nu)*m.cos(phi)**2
    denum=m.cos(phi)**2-m.cos(i)**2

    ang = m.atan(num/m.sqrt(abs(denum))) #absolute value for when latitude is near maximum/minimum
    ang=180*ang/m.pi

    return np.abs(ang)


def extract_odl_astL1A(fn):
    f = open(fn, 'r')
    body = f.read()

    def get_odl_parenth_value(text_odl, obj_name):
        posobj = str.find(text_odl, obj_name)
        posval = str.find(text_odl[posobj + 1:len(text_odl)], 'VALUE')
        posparenthesis = str.find(text_odl[posobj +1 +posval:len(text_odl)], '(')
        posendval = str.find(text_odl[posobj +1+ posval + posparenthesis:len(text_odl)], ')')

        val = text_odl[posobj+posval+posparenthesis+2:posobj+posval+posparenthesis+posendval +1]

        return val

    def get_odl_quot_value(text_odl, obj_name):
        posobj = str.find(text_odl, obj_name)
        posval = str.find(text_odl[posobj + 1:len(text_odl)], 'VALUE')
        posquote =  str.find(text_odl[posobj + 1 + posval:len(text_odl)], '"')
        posendval = str.find(text_odl[posobj + posval + posquote + 2:len(text_odl)], '"')

        val = text_odl[posobj + posval +posquote+ 2:posobj + posval + +posquote + posendval +2]

        return val

    # get latitude
    lat_val = get_odl_parenth_value(body, 'GRingPointLatitude')
    lat_tuple = [float(lat_val.split(',')[0]), float(lat_val.split(',')[1]), float(lat_val.split(',')[2]), float(lat_val.split(',')[3])]

    # get longitude
    lon_val = get_odl_parenth_value(body, 'GRingPointLongitude')
    lon_tuple = [float(lon_val.split(',')[0]), float(lon_val.split(',')[1]), float(lon_val.split(',')[2]), float(lon_val.split(',')[3])]

    # get calendar date + time of day
    caldat_val = get_odl_quot_value(body, 'CalendarDate')
    timeday_val = get_odl_quot_value(body, 'TimeofDay')
    caldat = dt(year=int(caldat_val.split('-')[0]), month=int(caldat_val.split('-')[1]), day=int(caldat_val.split('-')[2]),
                               hour=int(timeday_val.split(':')[0]), minute=int(timeday_val.split(':')[1]),
                               second=int(timeday_val.split(':')[2][0:2]), microsecond=int(timeday_val.split(':')[2][3:6]) * 1000)

    # get cloud cover
    cloudcov_val = get_odl_quot_value(body, 'SceneCloudCoverage')
    cloudcov_perc = int(cloudcov_val)

    # get flag if bands acquired or not: band 1,2,3N,3B,4,5,6,7,8,9,10,11,12,13,14
    list_band = []
    band_attr = get_odl_quot_value(body, 'Band3N_Available')
    band_avail = band_attr[0:3] == 'Yes'
    list_band.append(band_avail)
    band_attr = get_odl_quot_value(body, 'Band3B_Available')
    band_avail = band_attr[0:3] == 'Yes'
    list_band.append(band_avail)

    range_band = range(1, 15)
    range_band.remove(3)
    for i in range_band:
        band_attr = get_odl_quot_value(body, 'Band' + str(i) + '_Available')
        band_avail = band_attr[0:3] == 'Yes'
        list_band.append(band_avail)

    band_tags = pd.DataFrame(data=list_band,
                             index=['band_3N', 'band_3B', 'band_1', 'band_2', 'band_4', 'band_5', 'band_6', 'band_7',
                                    'band_8', 'band_9', 'band_10', 'band_11', 'band_12', 'band_13', 'band_14'])

    # get scene orientation angle
    orient_attr = get_odl_quot_value(body, 'ASTERSceneOrientationAngle')
    # orient_angl = float(orient_attr)
    #some .met metadata files are incomplete for angles
    orient_angl = float(15.)

    return lat_tuple, lon_tuple, caldat, cloudcov_perc, band_tags, orient_angl
