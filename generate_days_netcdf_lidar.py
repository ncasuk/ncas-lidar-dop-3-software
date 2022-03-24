#!/usr/local/anaconda3/envs/chil_3_8/bin/python
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
mpl.use('Agg')  #Needs to be before import halo_python as pylab gets called in that
import halo_python
#import halo_python_v3

import datetime
from datetime import date
from pylab import date2num

print("Python version in code = ", sys.version)

# -------------------------
# Explain function usage
# Needs to be before main()
# -------------------------
def usage():

    print("\n")
    print("Command line structure")
    print("./generate_days_netcdf_lidar.py -s yyyymmdd (of start day) -e yyyymmdd (of finish day) -l lidar choice (1 = CL51 ceilometer, 2 = Halo")
    print("-c = calibration factor to apply -r = for Halo only, read raw (= 1) or cleaned data (= 2)")
    print("For CL51, calibration factor = 0 means that no calibrated files are produced")
    print ("For Halo, calibration factor is only applied to cleaned data. Set = 1 to apply no further calibration")


# ------------------------------------------------------------------------
#def main():
# ------------------------------------------------------------------------
if len(sys.argv) <= 1:  #no arguments provided
    usage()
    sys.exit(2)
try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:e:l:c:r:")
except getopt.GetoptError as err:	#Changed since v2
    # print help information and exit:
    print(str(err)) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)


# ------------------
# Default cal_factor 
# ------------------
cal_factor = 1.0 

# ------------------
# Default lidar_type  
# ------------------
lidar_type = 0	#If a lidar type isn't entered neither option will be called 

# ----------------------
# Default Halo file type  
# ----------------------
infile_type = 1	#If a Halo input file type isn't entered, only unprocessed data will be read  

# -------------------------
# Set default date as today
# -------------------------
today = date.today()
year  = today.year
month = today.month
day   = today.day


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
    elif o in ("-c", "--cal_factor"):
        cal_factor = float(a)
    elif o in ("-r", "--infile_type"):	#Only for Halo
        infile_type = int(a)

    else:
        assert False, "unhandled option"

print("Start day: ",startday)
print("End day: ",endday)
print("Lidar type: ",lidar_type) 
print("Cal. factor: ",cal_factor) 

startnum = int(date2num(datetime.datetime(int(startday[0:4]),int(startday[4:6]),int(startday[6:8]),0,0,0)))
endnum = int(date2num(datetime.datetime(int(endday[0:4]),int(endday[4:6]),int(endday[6:8]),0,0,0)))
print('Start, end Julian date = ',startnum,endnum)

if lidar_type == 1:
    for nday in range(startnum,endnum+1):
        datevals = halo_python.generate_netcdf_common(nday)
        print('\n')
        print('Date being generated = ',datevals[3])
        print('\n')

        halo_python.generate_netcdf_cl51(nday, cal_factor)

elif lidar_type ==2:
    for nday in range(startnum,endnum+1):
        datevals = halo_python.generate_netcdf_common(nday)
        print('\n')
        print('Date being generated = ',datevals[3])
        print('\n')

        halo_python.generate_netcdf_halo(nday, infile_type, cal_factor)

else:
    print("Lidar type should be 1 for CL51 or 2 for Halo")


