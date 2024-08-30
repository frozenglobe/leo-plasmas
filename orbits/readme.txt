IMPORT FILES
intervals.pkl - < 1000 km intervals, table with datetime (timezone aware), unix, julian date (generated in '1000-int.ipynb')
whi_int.pkl - intervals in which WHISPER data is available in datetime (timezone aware)
utc_jday_conv.py - file for converting datetime strings (yy-mm-ddThh:mm:ssZ) to unix time and julian dates

26464tle.txt - C2 TLE between 2010-10-01 and 2011-12-01 (I think)
*c22011.txt - position and velocity information of C2 manually extracted from CSA .cef between 2011-01-01 00:00:00 and 2011-01-03 08:05:00
    used in 'sgp4-orbits.ipynb' and 'testing-coordinates.ipynb'

dm-intervals-c2-240613-093916.csv - .csv file containing intervals from CLUSTER data mining


JUPYTER NOTEBOOKS
sgp4-orbits.ipynb - learning to propagate orbits with SGP4 and comparing with CSA definitive position data; plotting orbits
1000-int.ipynb - created dataframe of time intervals when C2 was less than 1000 km above the surface of the Earth. GSE and TEME positional data of C2 between 2010-10-24 and 2011-12-14. Plotted something probably for a sanity check.
debris-speeds.ipynb - investigating how fast debris actually moves! 


OUTPUT FILES
c2orb-11.png - coincident orbits (yay) propagated with SGP4 and CSA aux definitive position data