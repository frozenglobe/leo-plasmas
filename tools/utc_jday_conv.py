#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import re
from datetime import datetime, timezone
from sgp4.api import jday

# removes the T and Z between the day and the time and splits each of yy/mm/dd hh:mm:ss into (six) columns
def time_split(file):
    ts = re.split("[TZ,:-]+", file[0])[:-1]
    for i in range(1,file.shape[0]):
        ts = np.vstack([ts,re.split("[TZ,:-]+", file[i])[:-1]])
    return ts

# returns UTC timestamps for one column of dates & times (separated by T and Z)
def utc(file):
    ts = time_split(file)
    utc = np.zeros(file.shape[0])
    for i in range(file.shape[0]):
        ts_store = [int(x) for x in ts[i]]
        utc[i] = datetime(ts_store[0],ts_store[1],ts_store[2],ts_store[3],ts_store[4],ts_store[5],tzinfo=timezone.utc).timestamp()
    return utc

# returns UTC timestamps + difference for two columns of dates & times
def utc_diff(file):
    ts = time_split(file)
    utc = np.zeros((file.shape[0],3))
    for i in range(file.shape[0]):
        ts_store = [int(x) for x in ts[i]]
        utc[i,0] = datetime(ts_store[0],ts_store[1],ts_store[2],ts_store[3],ts_store[4],ts_store[5]).timestamp()
        utc[i,1] = datetime(ts_store[6],ts_store[7],ts_store[8],ts_store[9],ts_store[10],ts_store[11]).timestamp()
        utc[i,2] = (utc[i,1] - utc[i,0]) # difference
    return utc

# returns jday; set add = True if require jdate (sums fraction)
def jd(file,add):
    ts = time_split(file)
    if add == True: 
        jd = np.zeros(file.shape[0])
        for i in range(file.shape[0]):
            ts_store = [int(x) for x in ts[i]]
            jd[i] = np.sum(jday(ts_store[0],ts_store[1],ts_store[2],ts_store[3],ts_store[4],ts_store[5]),axis=0)
    else: 
        jd = np.zeros((file.shape[0],2))
        for i in range(file.shape[0]):
            ts_store = [int(x) for x in ts[i]]
            jd[i,:] = jday(ts_store[0],ts_store[1],ts_store[2],ts_store[3],ts_store[4],ts_store[5])
    return jd

