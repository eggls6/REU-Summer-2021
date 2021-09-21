#!/usr/bin/env python
# coding: utf-8


import pandas as pd               
import numpy as np
import math

import spiceypy as sp
import astropy.coordinates
import re
import sgp4.api as sg
import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
import os
import sys
from timeit import default_timer as timer
from astropy.time import Time
from astropy.time import TimeDelta
import datetime as dt


# #### Let us start with reading in the first 10000 LSST observations for Main Belt Asteroids (S1 in the Synthetic Solar System model)
path = "/data/projects/lsst/baseline_fbs_v1p7p1/"
dir_list = os.listdir(path)

dflist=[]

for d in dir_list:
    if(d[0:2]=='S1'):
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0000000.h5'))
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0010000.h5'))

# every dataframe looks like this
dflist[0]

# concatenate them into a single dataframe
dfin=pd.concat(dflist)

# we can sort the resulting dataframe by by FieldID
dfin.sort_values(['FieldID'], inplace=True)

# then grouping and counting is a little faster. It seems that FieldID 0 has the most observations
dfin.groupby(['FieldID']).count()['ObjID']


def icrf2radec(pos, deg=True):
    """Convert ICRF xyz to Right Ascension and Declination.
    Geometric states on unit sphere, no light travel time/aberration correction.
    
    Parameters:
    -----------
    pos ... real, dim=[n, 3], 3D vector of unit length (ICRF)
    deg ... True: angles in degrees, False: angles in radians
    Returns:
    --------
    ra ... Right Ascension [deg]
    dec ... Declination [deg]
    """
    norm=np.linalg.norm
    array=np.array
    arctan2=np.arctan2
    arcsin=np.arcsin
    rad2deg=np.rad2deg
    modulo=np.mod
    pix2=2.*np.pi
    
    if(pos.ndim>1):
        r=norm(pos,axis=1)
        xu=pos[:,0]/r
        yu=pos[:,1]/r
        zu=pos[:,2]/r
    else:
        r=norm(pos)
        xu=pos[0]/r
        yu=pos[1]/r
        zu=pos[2]/r
    
    phi=arctan2(yu,xu)
    delta=arcsin(zu)
    
    if(deg):
        ra = modulo(rad2deg(phi)+360,360)
        dec = rad2deg(delta)
    else:
        ra = modulo(phi+pix2,pix2)
        dec = delta
    
    return ra, dec


def radec2icrf(ra, dec, deg=True):
    """Convert Right Ascension and Declination to ICRF xyz unit vector.
    Geometric states on unit sphere, no light travel time/aberration correction.
    
    Parameters:
    -----------
    ra ... Right Ascension [deg]
    dec ... Declination [deg]
    deg ... True: angles in degrees, False: angles in radians
    Returns:
    --------
    x,y,z ... 3D vector of unit length (ICRF)
    """
    deg2rad=np.deg2rad
    array=np.array
    cos=np.cos
    sin=np.sin
    
    if(deg):
        a = deg2rad(ra)
        d = deg2rad(dec)
    else:
        a = array(ra)
        d = array(dec)
    
    cosd = cos(d)
    x = cosd*cos(a)
    y = cosd*sin(a)
    z = sin(d)
    
    return array([x, y, z])


field_ids = dfin['FieldID'].unique()

dates = dfin['FieldMJD'].unique()

print("Earliest field date: " + str(min(dates)) + " or 10-01-2022 at 23:39:19.663 UTC")
print("Last field date: " + str(max(dates)) + " or 10-30-2022 at 05:03:20.104 UTC")

with open('starlink tle.txt') as f:
    starlinks = f.read().splitlines() 

chunks = [starlinks[n:n + 3] for n in range(0, len(starlinks), 3)]

print("There are " + str(len(chunks)) + " satellites")


import warnings


def bearing(lat1, long1, lat2, long2):
    '''
    Finds the bearing from one (latitude, longitude) point to another
    
    Parameters:
    -----------
    lat1, long1 ... real [rad]: Initial point
    lat2, long2 ... real [rad]: End point
    
    Returns:
    --------
    bear ... real [rad]: Bearing    
    '''
    
    b2 = np.sin(long2-long1) * np.cos(lat2)
    b1 = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(long2-long1)
    bear = np.arctan2(b2, b1)
    
    return bear
 

def cross_track_distance(lat1, long1, lat2, long2, lat3, long3):
    '''
    Finds the cross track distance from a point to line
    
    Parameters:
    -----------
    lat1, long1 ... real [rad]: Initial point of line
    lat2, long2 ... real [rad]: End point of line
    lat3, long3 ... real [rad]: Third point
    
    Returns:
    --------
    dist ... real [rad]: Cross track distance 
    '''
    c1 = SkyCoord(ra=long1*u.radian, dec=lat1*u.radian, frame='icrs')  # First point
    c2 = SkyCoord(ra=long3*u.radian, dec=lat3*u.radian, frame='icrs')  # Third point

    sep = c1.separation(c2)
    
    bear_1_3 = bearing(lat1, long1, lat3, long3)
    bear_1_2 = bearing(lat1, long1, lat2, long2)

    dist = np.arcsin(np.sin(sep.radian) * np.sin((bear_1_3) - bear_1_2))
    
    return dist


warnings.filterwarnings('ignore')

for ch in chunks:
    sat_num = ch[0].strip()  # STARLINK-xx
    s = ch[1]
    t = ch[2]
    sat = sg.Satrec.twoline2rv(s, t)  # check each satellite against each field
    
    flag_2 = False

    if flag_2:
        continue
                
    for f in field_ids:
        ra_vals = []  # for satellites
        dec_vals = []  # for satellites
        sat_times = []
        
        data = dfin[dfin['FieldID']==f]  # all objects in the field

        ra_f = data['AstRA(deg)']  # RA of all objects
        dec_f = data['AstDec(deg)']  # Dec of all objects
                
        min_ra, max_ra = min(ra_f), max(ra_f)
        min_dec, max_dec = min(dec_f), max(dec_f)

        mjd_t = data['FieldMJD'].unique()[0]  # field time
        
        t = Time(mjd_t, format='mjd')  # convert to Time format
        t_change = TimeDelta(2.0, format='sec')
        
        flag = False
        
        count = 1
        
        # Propagate every 2.0 sec starting 30 sec before field time until 30 sec after
        t_before = t - TimeDelta(30.0, format='sec')
        
        while count <= 61:
            jd_t = t_before + 2400000.5  # convert mjd to jd for sgp4 calculation
            fr, whole = math.modf(float(str(jd_t)))  # fr = stuff after decimal
    
            e, r, v = sat.sgp4(float(str(jd_t)), round(fr, 12))
            x, y, z = r[0], r[1], r[2]
            length = np.sqrt(x**2 + y**2 + z**2)
            norm_coords = np.array([x/length, y/length, z/length])  # normalizing onto unit sphere
            ra, dec = icrf2radec(norm_coords)  # convert to ra and dec
            
            if ra != ra:  # if nan (sat can't be propagated)
                flag = True
                break  # don't bother looking at all the times in while loop

            if (min_ra - 1.0 <= ra <= max_ra + 1.0) and (min_dec - 1.0 <= dec <= max_dec + 1.0):
                ra_vals.append(ra)
                dec_vals.append(dec)
                sat_times.append(t_before.value)
    
            t_before += t_change
            count += 1
            
        if flag:
            print()
            print("Could not propagate " + sat_num)
            print()
            flag_2 = True
            break  # break out of for loop of the fields, go to next satellite
                
        # only successful satellite points
        if any(ra_vals):
            good_obj_ra, good_obj_dec = [], []
            
            obj_count = 0
            
            # for satellites
            min_ra, max_ra = min(ra_vals), max(ra_vals)
            min_dec, max_dec = min(dec_vals), max(dec_vals)
            
            min_long, max_long = np.deg2rad(min_ra), np.deg2rad(max_ra)
        
            min_lat, max_lat = np.deg2rad(min_dec), np.deg2rad(max_dec)        

            for obj_r, obj_d in zip(ra_f, dec_f):
                obj_long = np.deg2rad(obj_r)
                obj_lat = np.deg2rad(obj_d)

                dist = cross_track_distance(min_lat, min_long, max_lat, max_long, obj_lat, obj_long)
                bound_d = 2.0 * np.pi/(180*3600)  # 2 arcseconds in radians
                
                if abs(dist) <= bound_d:
                    good_obj_ra.append(obj_r)
                    good_obj_dec.append(obj_d)
                    obj_count += 1
                    
            # only successful object points 
            if any(good_obj_ra):
                plt.figure(dpi=150,figsize=(4,4))
                plt.scatter(ra_vals, dec_vals, s=1.5, label='Satellite')
                plt.scatter(good_obj_ra, good_obj_dec, s=1.5, label='Objects')
                title = "Field " + str(f) + " " + sat_num
                plt.title(title)
                plt.xlabel('RA (deg)')
                plt.ylabel('Dec (deg)')
                plt.legend()
#                 plt.savefig(title)
                plt.show()

                print("Number of objects affected: " + str(obj_count))
                print()
                print("Satellite times:")
                print(str(sat_times))
#                 print(str(sat_times)[1:-1])  # remove brackets from list output
