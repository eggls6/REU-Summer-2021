#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd               
import numpy as np
import math


# In[2]:


import spiceypy as sp
import astropy.coordinates
import re
import sgp4.api as sg
import astropy.units as u
from astropy.coordinates import SkyCoord


# In[3]:


import matplotlib.pyplot as plt
import os
import sys
from timeit import default_timer as timer
from astropy.time import Time
from astropy.time import TimeDelta
import datetime as dt


# #### Let us start with reading in the first 10000 LSST observations for Main Belt Asteroids (S1 in the Synthetic Solar System model)

# In[4]:


path = "/data/projects/lsst/baseline_fbs_v1p7p1/"
dir_list = os.listdir(path)

dflist=[]

for d in dir_list:
    if(d[0:2]=='S1'):
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0000000.h5'))
        dflist.append(pd.read_hdf('/data/projects/lsst/baseline_fbs_v1p7p1/'+d+'/visit-0010000.h5'))


# In[5]:


# every dataframe looks like this
dflist[0]


# In[6]:


# concatenate them into a single dataframe
dfin=pd.concat(dflist)


# In[7]:


# we can sort the resulting dataframe by by FieldID
dfin.sort_values(['FieldID'], inplace=True)


# In[8]:


# then grouping and counting is a little faster. It seems that FieldID 0 has the most observations
dfin.groupby(['FieldID']).count()['ObjID']


# In[9]:


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


# In[10]:


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


# In[11]:


field_ids = dfin['FieldID'].unique()


# In[12]:


dates = dfin['FieldMJD'].unique()


# In[13]:


print("Earliest field date: " + str(min(dates)) + " or 10-01-2022 at 23:39:19.663 UTC")
print("Last field date: " + str(max(dates)) + " or 10-30-2022 at 05:03:20.104 UTC")


# In[14]:


# open the TLE files
with open('./starlink_tle.txt') as f:
    starlinks = f.read().splitlines() 


# In[15]:


# split up the list into 3-line lists for each satellite
chunks = [starlinks[n:n + 3] for n in range(0, len(starlinks), 3)]


# In[16]:


print("There are " + str(len(chunks)) + " satellites")


# In[17]:


import warnings


# In[18]:


def crosstrackdistance(lat1, long1, lat2, long2, lat3, long3):
    c1 = SkyCoord(ra=long1*u.radian, dec=lat1*u.radian, frame='icrs')  
    c2 = SkyCoord(ra=long2*u.radian, dec=lat2*u.radian, frame='icrs') 
    c3 = SkyCoord(ra=long3*u.radian, dec=lat3*u.radian, frame='icrs')
    
    P1 = radec2icrf(c1.ra.radian, c1.dec.radian, deg=False)
    P2 = radec2icrf(c2.ra.radian, c2.dec.radian, deg=False) 
    P3 = radec2icrf(c3.ra.radian, c3.dec.radian, deg=False) 
    
    S = np.cross(P1,P2)
    dot = np.dot(S,P3)
    
    dist = np.arcsin(dot/np.linalg.norm(S)) 
    
    return dist


# In[19]:


# Time calculations output warnings --> ignore them
warnings.filterwarnings('ignore')

# first iterate through each satellite
for ch in chunks[:7]:
    sat_num = ch[0].strip()  # the label for each satellite --> STARLINK-xx
    s = ch[1]  # first line of TLE set
    t = ch[2]  # second line of TLE set
    sat = sg.Satrec.twoline2rv(s, t)  
    
    flag_2 = False

    # if True, satellite can't be propagated, so continue to next satellite
    if flag_2:
        continue
    
    # check each satellite against each field
    for f in field_ids[100:1000]:
        ra_vals = []  # for satellites' right ascensions
        dec_vals = []  # for satellites' declinations
        sat_times = []  # for satellites' times
        
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
        
        # Start 30 seconds before field time
        t_before = t - TimeDelta(30.0, format='sec')
        
        # We're going to propagate every 2 seconds from 30 seconds before to 30 seconds after field time because exposure for the field is 30 seconds
        # Because of this, there will be 61 total iterations of propagations
        # If a satellite is within this time, it will leave a streak in the field image
        while count <= 61:
            jd_t = t_before + 2400000.5  # convert MJD to JD for sgp4 propagation calculation
            fr, whole = math.modf(float(str(jd_t)))  # fr = digits after decimal of MJD
    
            e, r, v = sat.sgp4(float(str(jd_t)), round(fr, 12))  # r is [x,y,z] for propagated satellite
            x, y, z = r[0], r[1], r[2]
            length = np.sqrt(x**2 + y**2 + z**2)
            norm_coords = np.array([x/length, y/length, z/length])  # normalizing onto unit sphere
            ra, dec = icrf2radec(norm_coords)  # convert to RA and Dec
            
            if ra != ra:  # if nan --> satellite can't be propagated to that time
                flag = True
                break  # don't bother looking at all the times --> break while loop (REDUNDANT)

            # if the object is in the field i.e. within 1 degree
            if (min_ra - 1.0 <= ra <= max_ra + 1.0) and (min_dec - 1.0 <= dec <= max_dec + 1.0):
                ra_vals.append(ra)
                dec_vals.append(dec)
                sat_times.append(t_before.value)
    
            t_before += t_change  # Propagate 2.0 seconds
            count += 1 
            
        if flag:
            print()
            print("Could not propagate " + sat_num)
            print()
            flag_2 = True
            break  # break out of for loop of the fields
                
        # Only satellite points that are within the field of objects
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

                bound_d = 2.0 * np.pi/(180*3600)  # 2 arcseconds in radians

                # calculate cross track distance between satellite streak and objects
                dist = crosstrackdistance(min_lat, min_long, max_lat, max_long, obj_lat, obj_long)
                
                # an object is considered obstructed if the satellite streak is within 2 arcseconds of it
                if abs(dist) <= bound_d:
                    good_obj_ra.append(obj_r)
                    good_obj_dec.append(obj_d)
                    obj_count += 1
                    
            # only plot objects that are obstructed 
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


# ### As seen in plot "Field 230 STARLINK-43", an object is being categorized as affected, even though visually it looks much farther away from the satellite streak than 2 arcseconds. This probable error is repeated several times throughout the execution of the code. It is not yet known if this is indeed an error in calculation, or in fact those objects are 2 arcseconds away from the satellite streak.

# In[ ]:




