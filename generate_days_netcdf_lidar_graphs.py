#!/usr/local/anaconda3/envs/chil_3_8/bin/python
#/usr/local/anaconda3/bin/python3
#/usr/bin/python
#Change this back to old v2.7 command if needed

# =========================================================================
# GENERATE DAYS NETCDF CL51 or Halo 
# Produces NetCDF data from a series of Vaisala CL51 ceilometer or Halo lidarprocessed files 
#
# Author:        Judith Jeffery, RAL 
# History:       Based on Python code written for Leosphere lidar by Chris Walden 
# Version:	 0.1
# Last modified: 09/08/17 
# =========================================================================

import getopt, sys
import matplotlib as mpl
mpl.use('Agg')

import lidar_plot

import datetime
from datetime import date
from pylab import date2num

# -------------------------
# Explain function usage
# Needs to be before main()
# -------------------------
def usage():

    print("\n")
    print("Command line structure")
    print("./generate_days_netcdf_lidar_graphs.py -s yyyymmdd (of start day) -e yyyymmdd (of finish day) -l lidar choice (1 = CL51 ceilometer, 2 = Halo)")
    print("-d raw or calibrated, (1 = raw data, 2 = calibrated data)")
    print("\n")


# ------------------------------------------------------------------------
#def main():
# ------------------------------------------------------------------------

if len(sys.argv) <= 1:  #no arguments provided
    usage()
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:e:l:d:")
except getopt.GetoptError as err:	#Changed since v2
    # print help information and exit:
    print(str(err)) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)


# ------------------
# Default lidar_type  
# ------------------
lidar_type = 1	#If a type isn't entered CL51 will be assumed 

# ------------------
# Default data_type  
# ------------------
data_type = 1	#If a type isn't entered raw data will be assumed 

# -------------------------
# Set default date as today
# -------------------------
today = date.today()
year  = today.year
month = today.month
day   = today.day
epoch_offset = 719163

# --------------------
# Command line parsing
# --------------------
for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-s", "--startday"):
        startday = str(a)
    elif o in ("-e", "--endday"):
        endday = str(a)
    elif o in ("-l", "--lidar"):
        lidar_type = int(a)
    elif o in ("-d", "--data"):
        data_type = int(a)

    else:
        assert False, "unhandled option"

print("Start day: ",startday)
print("End day: ",endday)
print("Lidar type: ",lidar_type) 

startnum = int(date2num(datetime.datetime(int(startday[0:4]),int(startday[4:6]),int(startday[6:8]),0,0,0)))
endnum = int(date2num(datetime.datetime(int(endday[0:4]),int(endday[4:6]),int(endday[6:8]),0,0,0)))
print('Start, end Julian date = ',startnum,endnum)

#Depending on version of Python, date2num can either give day wrt 01/01/0000 (earlier versions) or 01/01/1970 (later versions)
#This causes startnum, endnum to either be ~738000 or ~19000
#Other parts of code are based on 01/01/0000 so to accommodate different versions of Python, add to value if it's less than 700000

#if startnum < 100000:
    #startnum = startnum + epoch_offset
    #endnum = endnum + epoch_offset

if lidar_type == 1 or lidar_type == 2:
    for nday in range(startnum,endnum+1):
        if lidar_type == 1:
            lidar_plot.lidar_plotting(nday, lidar_type, data_type)
            data_type_cbh = 3
            lidar_plot.lidar_plotting_time_series(nday, lidar_type, data_type_cbh)
            data_type_cov = 4
            lidar_plot.lidar_plotting_time_series(nday, lidar_type, data_type_cov)
        else:
            lidar_plot.lidar_plotting(nday, lidar_type, data_type)

else:
    print("Lidar type should be 1 for CL51 or 2 for Halo")


