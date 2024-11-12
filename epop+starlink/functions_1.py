#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import re
import os
import time
import json
import requests
import configparser
from datetime import datetime, timedelta
from sgp4.api import Satrec, jday, days2mdhms
from kamodo_ccmc.flythrough.utils import ConvertCoord
from scipy.optimize import curve_fit
import cdflib

r_e = 6371.

# takes in datetimes or strings of the form 'YYYY-mm-dd hh:mm:ss'
def julian(time):
    j_date = list(map(int, str(time).split()[0].split('-')))
    j_time = list(map(int, str(time).split()[1].split(':')))
    j_dt = np.sum(jday(j_date[0], j_date[1], j_date[2], j_time[0], j_time[1], j_time[2]))
    
    return j_dt, j_date, j_time

# load CDF data
def cdf_info(file_path:str, fill_val:bool):
    file_in = cdflib.CDF(file_path)
    zvarinfo = file_in.cdf_info().zVariables
    
    znum = []
    zvar = []
    zdim = []
    zfill = []
    zdtype = []
    zunits = []

    for i in zvarinfo:
        ivar = file_in.varinq(i)
        znum.append(ivar.Num)
        zvar.append(ivar.Variable)
        zdim.append(ivar.Dim_Sizes)
        zdtype.append(ivar.Data_Type_Description)
        try:
            zfill.append(file_in.varattsget(i)['FILLVAL'])
            zunits.append(file_in.varattsget(i)['UNITS'])
        except KeyError:
            zfill.append('NaN')
            zunits.append('NaN')
            continue
    if (fill_val == True) & (len(fillval) == len(zvar)):
        return pd.DataFrame({'zvar': zvar, 'zdim': zdim, 'ztype': zdtype, 'zunits': zunits, 'zfv': zfill})
    else:
        return pd.DataFrame({'zvar': zvar, 'zdim': zdim, 'ztype': zdtype, 'zunits': zunits})


def clwhi_get(cl_id:str, start_date:str, end_date:str):
    start_date_rf = re.sub('-', '', start_date)
    end_date_rf = re.sub('-', '', end_date)
    base_fn = 'https://spdf.gsfc.nasa.gov/pub/data/cluster/' + cl_id + '/whi/efield_spectralpowerdensity_naturalmode/' + start_date_rf[:4] + '/' # obtain file names
    file_dir = 'data_psd' # creates new directory in current directory to store data downloads

    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
        print('Created', file_dir, 'directory')
    else:
        print('Directory', file_dir, 'already exists')

    a = str(requests.get(base_fn).content)
    b = re.sub(r'<.*?>', '', a)
    c = re.sub(r"b.*?Directory&nbsp;  - ", '',b)
    d = re.sub(r"\\n\\n\\n'", '', c)
    e = re.sub(r"20..-.*?M", '',d)
    f = e.split('\\n')[1:]
    
    index_1 = [idx for idx, s in enumerate(f) if start_date_rf in s][0]
    index_2 = [idx for idx, s in enumerate(f) if end_date_rf in s][0]
    
    for i in f[index_1:index_2+1]:
        file_path = os.path.join(file_dir, i)
        if os.path.exists(file_path):
            print('File', i, 'already exists')
            continue
        else:
            with open(file_path, 'wb') as file:
                print('Downloading', i)
                resp = requests.get(base_fn + i)
                file.write(resp.content)
    return

def clsta_get(cl_id:str, start_date:str, end_date:str):
    start_date_rf = re.sub('-', '', start_date)
    end_date_rf = re.sub('-', '', end_date)
    base_fn = 'https://spdf.gsfc.nasa.gov/pub/data/cluster/' + cl_id + '/sta/powerspectraldensity/' + start_date_rf[:4] + '/' # obtain file names
    file_dir = 'data_psd' # creates new directory in current directory to store data downloads

    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
        print('Created', file_dir, 'directory')
    else:
        print('Directory', file_dir, 'already exists')

    a = str(requests.get(base_fn).content)
    b = re.sub(r'<.*?>', '', a)
    c = re.sub(r"b.*?Directory&nbsp;  - ", '',b)
    d = re.sub(r"\\n\\n\\n'", '', c)
    e = re.sub(r"20..-.*?M", '',d)
    f = e.split('\\n')[1:]
    
    index_1 = [idx for idx, s in enumerate(f) if start_date_rf in s][0]
    index_2 = [idx for idx, s in enumerate(f) if end_date_rf in s][0]
    
    for i in f[index_1:index_2+1]:
        file_path = os.path.join(file_dir, i)
        if os.path.exists(file_path):
            print('File', i, 'already exists')
            continue
        else:
            with open(file_path, 'wb') as file:
                print('Downloading', i)
                resp = requests.get(base_fn + i)
                file.write(resp.content)
    return

# obtain data from MMS REST HTTPs API service. Dates in format YYYY-mm-dd. Example: descriptor: 'epsd' and data_rate: 'fast' 
def mms_get(mms_id:str, start_date:str, end_date:str, descriptor:str, data_rate:str):
    base_fn = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/file_names/science?' # obtain file names
    base_dl = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?' # download
    file_dir = 'data_' + descriptor # creates new directory in current directory to store data downloads
    query = 'sc_id=' + mms_id + '&start_date=' + start_date + '&end_date=' + end_date + '&descriptor=' + descriptor + '&data_rate_mode=' + data_rate
    
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
        print('Created', file_dir, 'directory')
    else:
        print('Directory', file_dir, 'already exists')

    a = requests.get(base_fn + query).content
    b = str(a)[2:-1].split(',')
    c = [i.split('/')[-1] for i in b]

    for i in c:
        file_path = os.path.join(file_dir, i)
        if os.path.exists(file_path):
            print('File', i, 'already exists')
            continue
        else:
            with open(file_path, 'wb') as file:
                print('Downloading', i)
                resp = requests.get(base_dl + 'file=' + i)
                file.write(resp.content)
    return

# t1 and t2 are dates in the form 'YYYY-mm-dd'
def findTLE(sat, t1, t2, config_path:str):
    config = configparser.ConfigParser()
    config.read(config_path)
    configUsr = config.get("configuration","username")
    configPwd = config.get("configuration","password")
    siteCred = {'identity': configUsr, 'password': configPwd}
    requestFindTLE   = "/class/gp_history/NORAD_CAT_ID/" + str(sat) + "/EPOCH/" + str(t1) + "--" + str(t2) + "/predicates/EPOCH,TLE_LINE1,TLE_LINE2/format/json/"
    
    with requests.Session() as session:
        resp = session.post('https://www.space-track.org/ajaxauth/login', data = siteCred)
        if resp.status_code != 200:
            raise MyError(resp, "POST fail on login")

        resp = session.get('https://www.space-track.org/basicspacedata/query' + requestFindTLE)
        if resp.status_code != 200:
            print(resp)
            raise MyError(resp, "GET fail on request for satellites")
        
        epoch = []
        tle_1 = []
        tle_2 = []
        retData = json.loads(resp.text)
        for t in retData:
            epoch.append(t['EPOCH'])
            tle_1.append([t['TLE_LINE1']])
            tle_2.append([t['TLE_LINE2']])
        
    session.close()
    
    a = [datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f') for i in epoch]
    tle_df = pd.DataFrame({'tle1': tle_1, 'tle2': tle_2}, index=a)
    
    return tle_df

def tle_get(epoch_start, days_no, launch_dt:str, decay_dt:str, periapsis): 
    epoch_range = []
    for j in range(days_no + 1):
        epoch_range.append((epoch_start + timedelta(days=j)).strftime('%Y-%m-%d'))

    with requests.Session() as session:

        resp = session.post(uriBase + requestLogin, data = siteCred)
        if resp.status_code != 200:
            raise MyError(resp, "POST fail on login")

        satIds = []
        epoch = []
        tle_1 = []
        tle_2 = []
        
        count=1
        
        for i in range(len(epoch_range)-1):
            if count > 10:
                print("Snoozing for 90 secs for rate limit reasons (max 30/min and 300/hr)...")
                time.sleep(90)
                count = 1
            
            epoch_s = epoch_range[i]
            epoch_e = epoch_range[i+1]

            requestTLE_1 = "/class/gp_history/LAUNCH_DATE/<" + launch_dt + "/DECAY_DATE/null-val/OBJECT_TYPE/PAYLOAD/PERIAPSIS/<" + str(periapsis) + "/EPOCH/" + epoch_s + "--" + epoch_e + "/predicates/NORAD_CAT_ID,EPOCH,TLE_LINE1,TLE_LINE2/format/json"
            resp = session.get(uriBase + requestCmdAction + requestTLE_1)
            if resp.status_code != 200:
                print(resp)
                raise MyError(resp, "GET fail on request for satellites")
            
            retData = json.loads(resp.text)
            for e in retData:
                satIds.append(e['NORAD_CAT_ID'])
                tle_1.append([e['TLE_LINE1']])
                tle_2.append([e['TLE_LINE2']])
                epoch.append(e['EPOCH'])
            
            requestTLE_2 = "/class/gp_history/LAUNCH_DATE/<" + launch_dt + "/DECAY_DATE/>" + decay_dt + "/OBJECT_TYPE/PAYLOAD/PERIAPSIS/<" + str(periapsis) + "/EPOCH/" + epoch_s + "--" + epoch_e + "/predicates/NORAD_CAT_ID,EPOCH,TLE_LINE1,TLE_LINE2/format/json"
            resp = session.get(uriBase + requestCmdAction + requestTLE_1)
            if resp.status_code != 200:
                print(resp)
                raise MyError(resp, "GET fail on request for satellites")
            
            retData = json.loads(resp.text)
            for e in retData:
                satIds.append(e['NORAD_CAT_ID'])
                tle_1.append([e['TLE_LINE1']])
                tle_2.append([e['TLE_LINE2']])
                epoch.append(e['EPOCH'])
            
            count+=1
            
    a = [datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f') for i in epoch]
    tle_df = pd.DataFrame({'id': satIds, 'epoch': a, 'tle_1': tle_1, 'tle_2': tle_2})

    session.close()
    
    return tle_df

def find_conj(dt_start, days_no, obs, min_buffer, is_dist_buffer, peri_dt_list):
    tle_df = pd.read_pickle('tle_' + str(dt_start.year) + '_' + ('%02d' % dt_start.month) + '.pkl')
    
    is_mindist_list = []
    is_mindist_t = []
    is_mindist_id = []
    is_mindist_alt = []
    
    peri_dt_trun = peri_dt_list[peri_dt_list > dt_start][peri_dt_list < (dt_start + timedelta(days=days_no-1))]
    
    for i in peri_dt_trun:
        j_sta = julian(i - timedelta(minutes=min_buffer))[0]
        j_end = julian(i + timedelta(minutes=min_buffer))[0]
        j_ran = np.arange(j_sta, j_end, 1/86400)[:-1]
        
        d = [*set(tle_df['id'])] # all sat norads
        e = tle_df.set_index('id')
        
        # observation satellite
        obs_e = e.loc[obs]['epoch']
        obs_i = np.argmin([abs(j.total_seconds()) for j in (obs_e - i)])
        obs_s = e.loc[obs].iloc[obs_i]['tle_1'][0]
        obs_t = e.loc[obs].iloc[obs_i]['tle_2'][0]

        d.remove(obs)

        for j in d:
            f = e.loc[j]['epoch']
            g = np.argmin([abs(j.total_seconds()) for j in (f - i)])
            s = e.loc[j].iloc[g]['tle_1'][0]
            t = e.loc[j].iloc[g]['tle_2'][0]

            h = Satrec.twoline2rv(s, t).sgp4_array(j_ran, np.zeros(len(j_ran)))[1]
            k = Satrec.twoline2rv(obs_s, obs_t).sgp4_array(j_ran, np.zeros(len(j_ran)))[1]     
            is_dist = np.sum( (h - k)**2, axis=1 )**0.5

            is_mindist = np.min(is_dist)

            if is_mindist < is_dist_buffer:
                is_mindist_list.append(is_mindist)
                sec_diff = int(np.argmin(is_dist) - int(len(j_ran)/2))
                is_mindist_t.append(i + timedelta(seconds = sec_diff))
                is_mindist_id.append(j)
                alt = np.sum(k[np.argmin(is_dist)]**2)**0.5 - r_e
                is_mindist_alt.append(alt)
            else:
                continue
        
    return is_mindist_id, is_mindist_t, is_mindist_list, is_mindist_alt

def tle_select(ts_list, time):
    TLE_epoch = []
    for i in range(len(ts_list)):
        tle_time = ts_list['tle1'].iloc[i][0].split()[3]
        a = 2000 + int(tle_time[:2])
        b, c, d, e, f = days2mdhms(int(tle_time[:2]), float(tle_time[2:]))
        j_tle = np.sum(jday(a, b, c, d, e, f))
        TLE_epoch = np.append(TLE_epoch, j_tle)
        arg = np.argmin(abs(TLE_epoch - time))

    return ts_list.iloc[arg]

# B-field vector generation from IGRF model. Only takes in GSE.
def igrf(time, x, y, z, conj_time):
    lon_gdz, lat_gdz, alt_gdz = ConvertCoord(time, x/r_e, y/r_e, z/r_e, 'GSE', 'car', 'GDZ', 'sph', verbose=False)[:-1]
    
    # setting up lists
    mag = np.zeros(len(lon_gdz))
    diff_lat = np.zeros(len(lon_gdz))
    diff_lon = np.zeros(len(lon_gdz))
    diff_alt = np.zeros(len(lon_gdz))
    diff_alt_km = 0.001 # in km, i.e. 1 m
    eff_r = alt_gdz + r_e # effective radius
    
    # obtain IGRF field data and generate differentials
    for i in range(len(lon_gdz)):
        r = requests.get('https://geomag.bgs.ac.uk/web_service/GMModels/igrf/13/?latitude=' + str(lat_gdz[i]) + '&longitude=' + str(lon_gdz[i]) + '&altitude=' + str(alt_gdz[i]) + '&date=' + str(conj_time.date()) + '&format=json')
        load_field = json.loads(r.text)['geomagnetic-field-model-result']['field-value']
        vert = load_field['vertical-intensity']['value']
        n_prop = load_field['north-intensity']['value'] / abs(vert)
        e_prop = load_field['east-intensity']['value'] / abs(vert)
        diff_lat[i] = n_prop / (2 * np.pi * eff_r[i] / diff_alt_km) * 360
        diff_lon[i] = e_prop / (2 * np.pi * eff_r[i] * np.cos(lat_gdz[i]/180*np.pi) / diff_alt_km) * 360
        diff_alt[i] = vert / abs(vert) * diff_alt_km
        mag[i] = load_field['total-intensity']['value']
    
    # second position to convert back into GSE
    sec_lat = lat_gdz + diff_lat
    sec_lon = lon_gdz + diff_lon
    
    x_2, y_2, z_2 = ConvertCoord(time, sec_lon, sec_lat, alt_gdz - diff_alt, 'GDZ', 'sph', 'GSE', 'car')[:-1]
    u_gse = x_2*r_e - x
    v_gse = y_2*r_e - y
    w_gse = z_2*r_e - z
    
    field_vec = np.stack((u_gse, v_gse, w_gse), axis=1)
    
    return field_vec

# Quadratic trajectory fitting
def model(t, a, b, c):
    return a*t**2 + b*t + c

def curve_gen(time, pmatrix):
    return model(time, pmatrix[0], pmatrix[1], pmatrix[2])

def fit_gen(r, s, time, sat_pos, check:bool):
    
    # The time parameter requires some moulding in order to converge to the correct solution. 
    # For example, the unix times are typically too big and hence require a rounded subtraction. 
    # The number to subtract is completed by scaling the time down by 'r' and then forcing it to 
    # be an integer - this is what I refer to by 'rounding'. Then 's' is the scaling. The new
    # time variable is 't_fit'. 
    t_sub = int(time[0]*r)/r
    t_fit = (time - t_sub) * s

    px, ex = curve_fit(model, t_fit, sat_pos[:,0])
    py, ey = curve_fit(model, t_fit, sat_pos[:,1])
    pz, ez = curve_fit(model, t_fit, sat_pos[:,2])
    
    p_fit = np.array((px, py, pz)).T
        
    # Checking vaildity of curve fit. The first ratio should be much smaller than one, 
    # and the condition number should be as small as possible. 
    if check == True: 
        print('Ratio of covariance diagonal to parameter:', np.max(np.diag(ex)/px))
        print('Condition number of covariance matrix:', round(np.max([np.linalg.cond(ex), np.linalg.cond(ey), np.linalg.cond(ez)])))
    
    x_fit = curve_gen(t_fit, px)
    y_fit = curve_gen(t_fit, py)
    z_fit = curve_gen(t_fit, pz)
    
    r_fit = np.stack((x_fit, y_fit, z_fit), axis=1)
    
    return r_fit, p_fit, np.array([t_sub, s]), t_fit

# Surface-curve intersection search
def is_search(sat_normal, sat_pos, obs_pos, field_vec, buffer):
    find_list = []
    find_obs_index = []
    find_sat_index = []
    mu_dist_list = []
        
    for j in range(len(sat_normal)):
        find_list_n = []
        find_obs_index_n = []
        mu_dist_list_n = []
        field_vec_no = field_vec[j]
        sat_norm_no = sat_normal[j]
        
        for i in range(len(obs_pos)):
            int_vec_diff = obs_pos[i] - sat_pos[j] # vector of intersection - point
            cross = np.cross(field_vec_no, sat_norm_no)
            cross_norm = cross/np.sum(cross**2)**0.5
            mu_dist = abs(np.dot(int_vec_diff, cross_norm))
            
            if mu_dist < buffer:
                find_list_n.append(abs(np.dot(sat_norm_no, obs_pos[i]) - np.dot(sat_norm_no, sat_pos[j])))
                find_obs_index_n.append(i)
                mu_dist_list_n.append(mu_dist)
            else: 
                continue
        
        if len(find_list_n) > 0:
            find_list.append(np.min(find_list_n))
            find_obs_index.append(find_obs_index_n[np.argmin(find_list_n)])
            mu_dist_list.append(mu_dist_list_n[np.argmin(find_list_n)])
            find_sat_index.append(j)
    
    print(len(find_list), 'out of', len(sat_normal))
    
    vec_index = find_sat_index[np.argmin(find_list)]
    obs_pos_index = find_obs_index[np.argmin(find_list)]
    
    return vec_index, obs_pos_index, mu_dist_list, find_list