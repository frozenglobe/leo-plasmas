#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cdflib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sgp4.api import Satrec, jday, days2mdhms
from datetime import datetime
import re
r_e = 6371.

# returns the first TLE of the day
def TLEday(file):
    lines = open(file,'r').readlines()
    perday = lines[0:2]
    for i in range(1,int(len(lines)/2)): # compares the integer TLE epoch (the 4th field, 5 digits)
        if (lines[2*(i)].split()[3][:5] != lines[2*(i-1)].split()[3][:5]):
            perday = np.append(perday, [[lines[2*i]], [lines[2*i+1]]])
    return perday

# generates position array of the day of the TLE
def pos_day(tle1, tle2):
    d1 = 2000 + int(tle1.split()[3][:2])
    d2, d3, t1, t2, t3 = days2mdhms(int(tle1.split()[3][:2]),float(tle1.split()[3][2:]))
    jd_start = jday(d1, d2, d3, t1, t2, t3)[0]
    jd_end = jd_start+1
    jd_range = np.arange(jd_start, jd_end, 1/86400) # range of julian dates to generate, interval is one second
    sat = Satrec.twoline2rv(tle1,tle2) # reading TLE
    r = np.zeros([jd_range.shape[0],3]) # initialising table of positions
    for i in range(jd_range.shape[0]):
        r[i,:] = sat.sgp4(jd_range[i],0)[1] # obtaining table of x y z
    return r

# generates position array given start and end Julian dates, for a file of TLEs containing the interval of interest
def pos_file(file, jd_start, jd_end):
    perday = TLEday(file) # generating TLEs
        
    TLE_epoch = [] # identifying the correct TLE to use
    for i in range(int(len(perday)/2)):
        d1 = 2000 + int(perday[2*i].split()[3][:2])
        d2, d3, t1, t2, t3 = days2mdhms(int(perday[2*i].split()[3][:2]),float(perday[2*i].split()[3][2:]))
        julian = jday(d1, d2, d3, t1, t2, t3)
        TLE_epoch = np.append(TLE_epoch, julian[0])
    
    arg = np.argmin(abs(TLE_epoch-jd_start))

    jd_range = np.arange(jd_start, jd_end, 1/86400) # range of julian dates to generate, interval is one second
    sat = Satrec.twoline2rv(perday[2*arg],perday[2*arg+1]) # reading TLE
    r = np.zeros([jd_range.shape[0],3]) # initialising table of positions
    for i in range(jd_range.shape[0]):
        r[i,:] = sat.sgp4(jd_range[i],0)[1] # obtaining table of x y z
    return r

# generates position array given start and end Julian dates, for an array of TLEs containing the interval of interest
def pos_tle(tle, jd_start, jd_end):        
    TLE_epoch = [] # identifying the correct TLE to use
    for i in range(int(len(tle)/2)):
        d1 = 2000 + int(tle[2*i].split()[3][:2])
        d2, d3, t1, t2, t3 = days2mdhms(int(tle[2*i].split()[3][:2]),float(tle[2*i].split()[3][2:]))
        julian = jday(d1, d2, d3, t1, t2, t3)
        TLE_epoch = np.append(TLE_epoch, julian[0])
    
    arg = np.argmin(abs(TLE_epoch-jd_start))

    jd_range = np.arange(jd_start, jd_end, 1/86400) # range of julian dates to generate, interval is one second
    sat = Satrec.twoline2rv(tle[2*arg],tle[2*arg+1]) # reading TLE
    r = np.zeros([jd_range.shape[0],3]) # initialising table of positions
    for i in range(jd_range.shape[0]):
        r[i,:] = sat.sgp4(jd_range[i],0)[1] # obtaining table of x y z
    return r

# compares, for each timestep in each interval, the distance between the chosen satellite and C2 (26464); the number argument is used to specify how many to show
def conj_file(file, number):
    int_df = pd.read_pickle('intervals.pkl')
    min_dist = np.zeros(int_df.shape[0]) # initialises table of minimum distances in each interval
    min_dist_sec = np.zeros(int_df.shape[0]) # time of minimum distances in seconds after start of interval
    for i in range(int_df.shape[0]):
        abs_diff = abs(pos_file('26464tle.txt', int_df['dt1-jd'][i], int_df['dt2-jd'][i]) - pos_file(file, int_df['dt1-jd'][i], int_df['dt2-jd'][i]))
        min_dist[i] = np.min( (np.sum(abs_diff**2,axis=1) )**0.5)
        min_dist_sec[i] = np.argmin( (np.sum(abs_diff**2,axis=1) )**0.5)
    
    minsorted = np.zeros(number)
    minsecsorted = np.zeros(number)
    arg = np.zeros(number)
    for i in range(number):
        minsorted[i] = np.sort(min_dist)[i]
        arg[i] = np.argmin(abs(min_dist - minsorted[i]))
        minsecsorted[i] = min_dist_sec[int(arg[i])]
    
    conj = pd.DataFrame({'int-start': int_df['dt1'][arg], 'min-s': minsecsorted, 'min-km': minsorted})
    return conj

def conj_tle(tle, number):
    int_df = pd.read_pickle('intervals.pkl')
    min_dist = np.zeros(int_df.shape[0]) # initialises table of minimum distances in each interval
    min_dist_sec = np.zeros(int_df.shape[0]) # time of minimum distances in seconds after start of interval
    for i in range(int_df.shape[0]):
        abs_diff = abs(pos_file('26464tle.txt', int_df['dt1-jd'][i], int_df['dt2-jd'][i]) - pos_tle(tle, int_df['dt1-jd'][i], int_df['dt2-jd'][i]))
        min_dist[i] = np.min( (np.sum(abs_diff**2,axis=1) )**0.5)
        min_dist_sec[i] = np.argmin( (np.sum(abs_diff**2,axis=1) )**0.5)
    
    minsorted = np.zeros(number)
    minsecsorted = np.zeros(number)
    arg = np.zeros(number)
    
    for i in range(number):
        minsorted[i] = np.sort(min_dist)[i]
        arg[i] = np.argmin(abs(min_dist - minsorted[i]))
        minsecsorted[i] = min_dist_sec[int(arg[i])]
    
    conj = pd.DataFrame({'int-start': int_df['dt1'][arg], 'min-s': minsecsorted, 'min-km': minsorted})
    return conj

# ub = upper bound
def conj_tle_ub(tle, ub):
    int_df = pd.read_pickle('intervals.pkl')
    min_dist = np.zeros(int_df.shape[0]) # initialises table of minimum distances in each interval
    min_dist_sec = np.zeros(int_df.shape[0]) # time of minimum distances in seconds after start of interval
    for i in range(int_df.shape[0]):
        abs_diff = abs(pos_file('26464tle.txt', int_df['dt1-jd'][i], int_df['dt2-jd'][i]) - pos_tle(tle, int_df['dt1-jd'][i], int_df['dt2-jd'][i]))
        min_dist[i] = np.min( (np.sum(abs_diff**2,axis=1) )**0.5)
        min_dist_sec[i] = np.argmin( (np.sum(abs_diff**2,axis=1) )**0.5)
    
    number = len(np.nonzero(min_dist < ub)[0])
    minsorted = np.zeros(number)
    minsecsorted = np.zeros(number)
    arg = np.zeros(number)
     
    for i in range(number):
        minsorted[i] = np.sort(min_dist)[i]
        arg[i] = np.argmin(abs(min_dist - minsorted[i]))
        minsecsorted[i] = min_dist_sec[int(arg[i])]
    
    conj = pd.DataFrame({'int-start': int_df['dt1'][arg], 'min-s': minsecsorted, 'min-km': minsorted})
    return conj