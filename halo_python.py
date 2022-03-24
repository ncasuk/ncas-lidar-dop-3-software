
# =========================================================================
# HALO PYTHON MODULE
# Tools to produces NetCDF data from a day of Halo lidar processed files
#
# Author:        Judith Jeffery, RAL 
# History:       Based on Python code written for Leosphere lidar by Chris Walden 
# Version:	 1.0
# Last modified: 09/08/17
# =========================================================================
module_version = 1.0

# -------------------------------------------------------------------------
# Import required tools
# -------------------------------------------------------------------------
import numpy as np
import os, re, sys, getopt, shutil, zipfile, string, pwd
import netCDF4 as nc4

import datetime
import time 
import calendar
import scipy.signal
from pylab import *
import module_data_object_python3
import lidar_plot

# ----------------------
# Define useful function
# ----------------------
def in_interval(seq,xmin,xmax):
    for i, x in enumerate(seq):
        if x>=xmin and x<xmax:
            yield i

# ------------------------------------------------------------------------
# Define the main netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_halo(nday,infile_type,cal_factor):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------
    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]

    #Get previous day for checking for data in previous day's file
    prev_datevals=generate_netcdf_common((nday-1))
    prev_datestring = prev_datevals[3]

    dat              = datetime.datetime(year,month,day,0,0,0)
    start_of_day     = date2num(dat)
    daystr           = dat.strftime("%Y%m%d_")
    dirstr           = dat.strftime("%Y") + '/' + dat.strftime("%Y%m") + '/' + dat.strftime("%Y%m%d") + '/'
    start_of_day_str = dat.strftime("%Y-%m-%d %H:%M")
    dat_str          = dat.strftime("%Y-%m-%dT%H:%M:%S") 

    print("start_of_day, daystr, dirstr, datestring, start_of_day_str = ", start_of_day, daystr, dirstr, datestring, start_of_day_str) 

    # -------------------------
    # Define various parameters
    # -------------------------
    max_h      = 9.6                   # Max height (km)
    hbin       = 0.03                  # Bin height (km)
    nh         = 1000		       # Initial number of height bins, trimmed later
    nangles    = 1			# For index_of_angles
    nangles_vad = 24 
    missing_value = -1.0E+20
    uncertainty_threshold = 1.0 #Uncertainty in cleaned data above which qc_flag is set not equal to 1

    # ----------------------
    # Initialize data arrays
    # ----------------------
    vals_co     = np.zeros((50000,nh,3))
    vals_cr     = np.zeros((50000,nh,3))
    vals_vad     = np.zeros((50000,nh,3))
    timesecs_co   = np.zeros((50000))
    timesecs_cr   = np.zeros((50000))
    timesecs_vad   = np.zeros((50000))
    elev_angle_co = np.zeros((50000))
    elev_angle_cr = np.zeros((50000))
    elev_angle_vad = np.zeros((50000,nangles_vad))
    azi_angle_vad = np.zeros((50000,nangles_vad))
    ok_to_write = np.zeros((infile_type,3))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/lidar-doppler-halo/mirror_lidar-doppler-halo2_data/Proc/" + dirstr 
    path_in_cr = "/data/lidar-doppler-halo/mirror_lidar-doppler-halo2_data/Proc/cross/" + dirstr
    cleaned_nc_path_in = "/data/lidar-doppler-halo/netcdf/calibrated/halo-doppler-lidar-118-co/"
    cleaned_vad_nc_path_in = "/data/lidar-doppler-halo/netcdf/calibrated/halo-doppler-lidar-118-vad75/" + str(year) + '/' + datestring + '/'
    cfarr_head = "ncas-lidar-dop-3_cao_"
    path_out = "/data/amof-netCDF/ncas-lidar-dop-3/"
    text_path_out = "/mnt/wilma_home/jla/halo/text_profile/"
    errorfile = path_out + "halo_errors.txt"
    nproffile = path_out + "number_of_halo_profiles.txt"
    midnight_file_root = path_out + "halo_midnight_"
    oformaterr1="%s"
    oformaterr2="%2i"
    oformaterr3="%5i"
    oformatmid1="%13.6e"

    if os.access(path_in, os.F_OK):	#A data directory exists for the day
        files = os.listdir(path_in)
        #tempstr = "Stare_118_"+datestring+'_'
        #expr = re.compile(tempstr+'[0-2][0-9]\.hpl')
        expr = re.compile('hpl')
        infiles = [elem for elem in files if expr.search(elem)]
        infiles.sort()
        #print(infiles)
        nfiles = len(infiles)
        if nfiles == 0:
            print('No co-polar files to process')

        files_cr = os.listdir(path_in_cr)
        infiles_cr = [elem for elem in files_cr if expr.search(elem)]
        infiles_cr.sort()
        #print infiles_cr
        nfiles_cr = len(infiles_cr)
        if nfiles_cr == 0:
            print('No cross-polar files to process')
    else:
        nfiles = 0
        nfiles_cr = 0
        nfiles_vad = 0
        print('No co- or cross- polar files to process - no directory for that day')

    # ----------------------
    # Process co-polar files
    # ----------------------
    print('Starting to process co-polar files')
    n_co = 0	#Times
    m_co = 0	#Range gates
    n_vad = 0

    #Check for any data from the previous day
    midnight_file = midnight_file_root + prev_datestring + '_0.txt'
    #Check this file was created
    if os.path.isfile(midnight_file):
        mid_file = open(midnight_file,'r') 
        #Check date in 1st line of file. May need to say how many records are in file - assuming 1 for now 
        line = mid_file.readline()
        verify_data = line.split()
        #data_type = verify_data[1]
        read_ngates = int(verify_data[2]) 
        print('Getting co-polar data from previous day')
        timesecs_co[n_co] = float(verify_data[3])/3600.0
        elev_angle_co[n_co] = float(verify_data[4])/3600.0
        for z in range(read_ngates):	#Put value in midnight file
        #for z in range(320):	#Put value in midnight file
            line = mid_file.readline()
            fdata = line.split()
            if len(fdata) > 0:
                vals_co[n_co,m_co,0] = float(fdata[0])	#Beta
                vals_co[n_co,m_co,1] = float(fdata[1])	#Speed
                vals_co[n_co,m_co,2] = float(fdata[2])	#SNR+1
       	        m_co += 1

        m_co = 0	#Range bin counter, set back to zero 
        n_co += 1 

        mid_file.close()

    for nf in range(nfiles):	#Indent from here
        #Is this a Stare or VAD file?
        if np.char.find(infiles[nf],"Stare") >= 0:
            fixed = 1
        else:
            fixed = 0

        f = open(path_in+infiles[nf], 'r')

        for z in range(2):
            line = f.readline()
        #It would be nice to read all these constants into 1 array using a loop, but they have different types!
        line = f.readline()
        tempvar = line.split(':')
        ngates = int(tempvar[1])
        line = f.readline()
        tempvar = line.split(':')
        rgate = float(tempvar[1])
        line = f.readline()
        tempvar = line.split(':')
        gatelength = int(tempvar[1])
        line = f.readline()
        tempvar = line.split(':')
        pulsesray = int(tempvar[1])
        line = f.readline()	#No. rays in file
        line = f.readline()	#Scan type
        line = f.readline()	#Focus range
        tempvar = line.split(':')
        focusrange = int(tempvar[1])
        line = f.readline()
        line = f.readline()
        tempvar = line.split(':')
        res = float(tempvar[1])

        for z in range(6):	
            line = f.readline()

        # ---------------------------
        # Repeated data for each time
        # ---------------------------

        if fixed == 1:	#Stare file
            while True:

                line = f.readline()
                if not line: break
                m_co = 0	
                fdata = line.split()
                timesecs_co[n_co] = float(fdata[0])
                #elev_angle_co[n_co] = 90.0 - float(fdata[2])
                elev_angle_co[n_co] = float(fdata[2])
                #pitch_co[n_co] = float(fdata[3])
                #roll_co[n_co] = float(fdata[4])

                for z in range(ngates):
                    line = f.readline()
                    fdata = line.split()
                    if len(fdata) > 0:
                        vals_co[n_co,m_co,0] = float(fdata[3])	#Beta
                        vals_co[n_co,m_co,1] = float(fdata[1])	#Speed
                        vals_co[n_co,m_co,2] = float(fdata[2])	#SNR+1
       	                m_co += 1
                n_co += 1

        else:	#VAD file
            while True:

                line = f.readline()
                if not line: break
                m_co = 0	
                fdata = line.split()
                timesecs_vad[n_vad] = float(fdata[0])
                #elev_angle_co[n_co] = 90.0 - float(fdata[2])
                azi_angle_vad[n_vad,:] = float(fdata[1])
                elev_angle_vad[n_vad,:] = float(fdata[2])
                #pitch_co[n_co] = float(fdata[3])
                #roll_co[n_co] = float(fdata[4])

                for z in range(ngates):
                    line = f.readline()
                    fdata = line.split()
                    if len(fdata) > 0:
                        vals_vad[n_vad,m_co,0] = float(fdata[3])	#Beta
                        vals_vad[n_vad,m_co,1] = float(fdata[1])	#Speed
                        vals_vad[n_vad,m_co,2] = float(fdata[2])	#SNR+1
       	                m_co += 1
                n_vad += 1


        f.close()
    
    print('Number of co-polar, VAD profiles = ', n_co, n_vad)

    if n_co > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals_co = vals_co[0:n_co,0:ngates,:]	#I don't understand why this is 0:ngates and not 0:ngates-1
        timesecs_co   = timesecs_co[0:n_co]
        elev_angle_co = elev_angle_co[0:n_co]
        timesecs_co *= 3600.0

        #print('Co-polar times 0:9 = ', timesecs_co[0:9])
        #print('Elevation angles 0:9 = ', elev_angle_co[0:9])
        #print('Co-polar backscatter 0:4, 0:4 = ', vals_co[0:4,0:4,0])
        #print('Co-polar wind 0:4, 0:4 = ', vals_co[0:4,0:4,1])
        #print('Co-polar SNR+1 0:4, 0:4 = ', vals_co[0:4,0:4,2])

        nn_start = 0    #Co-polar
        new_ind = halo_check_monotonic(timesecs_co, vals_co, elev_angle_co, n_co, nn_start, datestring, ngates, errorfile, midnight_file_root)
        #print('co-polar new_ind = ',new_ind)
        if new_ind >= 0:        #Non-monotonic data found
            n_co = n_co - 1     #This assumes it's just the last point. Could perhaps make more general?
            timesecs_co = timesecs_co[0:new_ind]
            #azi_arr = azi_arr[0:new_ind,:]
            elev_angle_co = elev_angle_co[0:new_ind]
            vals_co = vals_co[0:new_ind,:,:]


    if n_vad > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals_vad = vals_vad[0:n_vad,0:ngates,:]	#I don't understand why this is 0:ngates and not 0:ngates-1
        timesecs_vad   = timesecs_vad[0:n_vad]
        azi_angle_vad = azi_angle_vad[0:n_vad,:]
        elev_angle_vad = elev_angle_vad[0:n_vad,:]
        timesecs_vad *= 3600.0

        #print('VAD time 0:20 = ', timesecs_vad[0:20])
        #print('VAD azimuth angle 0:20 = ', azi_angle_vad[0:20,0])
        #print('VAD elevation angle 0:20 = ', elev_angle_vad[0:9,0])
        #print('VAD backscatter 0:4, 0:4 = ', vals_vad[0:4,0:4,0])
        #print('VAD wind 0:4, 0:4 = ', vals_vad[0:4,0:4,1])
        #print('VAD SNR+1 0:4, 0:4 = ', vals_vad[0:4,0:4,2])

        nn_start = 2    #VAD
        new_ind = halo_check_monotonic(timesecs_vad, vals_vad, elev_angle_vad[:,0], n_vad, nn_start, datestring, ngates, errorfile, midnight_file_root)
        #print('VAD new_ind = ',new_ind)
        if new_ind >= 0:
            #print('Times at end of VAD file', timesecs_vad[(n_vad-24):n_vad])
            timesecs_vad[new_ind:n_vad] = timesecs_vad[new_ind:n_vad] + 86400.0
            #print('Corrected times at end of VAD file', timesecs_vad[(n_vad-24):n_vad])

    # -------------------------
    # Process cross-polar files
    # -------------------------

    print('Starting to process cross-polar files')
    n_cr = 0	#Times
    m_cr = 0	#Range gates

    #Check for any data from the previous day
    midnight_file = midnight_file_root + prev_datestring + '_1.txt'
    #Check this file was created
    if os.path.isfile(midnight_file):
        mid_file = open(midnight_file,'r') 
        #Check date in 1st line of file. May need to say how many records are in file - assuming 1 for now 
        line = mid_file.readline()
        verify_data = line.split()
        read_ngates = int(verify_data[2]) 
        print('Getting cross-polar data from previous day')
        timesecs_cr[n_cr] = float(verify_data[3])/3600.0
        elev_angle_cr[n_cr] = float(verify_data[4])/3600.0
        for z in range(read_ngates):	#Put value in midnight file
        #for z in range(320):	#Put value in midnight file
            line = mid_file.readline()
            fdata = line.split()
            if len(fdata) > 0:
                vals_cr[n_cr,m_cr,0] = float(fdata[0])	#Beta
                vals_cr[n_cr,m_cr,1] = float(fdata[1])	#Speed
                vals_cr[n_cr,m_cr,2] = float(fdata[2])	#SNR+1
       	        m_cr += 1

        m_cr = 0	#Range bin counter, set back to zero 
        n_cr += 1 

        mid_file.close()

    for nf in range(nfiles_cr):	#Indent from here
        f = open(path_in_cr+infiles_cr[nf], 'r')
        #Just read header lines this time, variables should be same as for co-polar files
        for z in range(17):
            line = f.readline()

        # ---------------------------
        # Repeated data for each time
        # ---------------------------

        while True:

            line = f.readline()
            if not line: break
            m_cr = 0	
            fdata = line.split()
            timesecs_cr[n_cr] = float(fdata[0])
            #elev_angle_cr[n_cr] = 90.0 - float(fdata[2])
            elev_angle_cr[n_cr] = float(fdata[2])
            #pitch_cr[n_cr] = float(fdata[3])
            #roll_cr[n_cr] = float(fdata[4])

            for z in range(ngates):
                line = f.readline()
                fdata = line.split()
                if len(fdata) > 0:
                    vals_cr[n_cr,m_cr,0] = float(fdata[3])	#Beta
                    vals_cr[n_cr,m_cr,1] = float(fdata[1])	#Speed
                    vals_cr[n_cr,m_cr,2] = float(fdata[2])	#SNR+1
                    m_cr += 1
            n_cr += 1

        f.close()


    # ---------------------
    # If data were read 
    # ---------------------
    if n_cr > 0:

        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals_cr = vals_cr[0:n_cr,0:ngates,:]
        timesecs_cr   = timesecs_cr[0:n_cr]
        elev_angle_cr = elev_angle_cr[0:n_cr]
        timesecs_cr   *= 3600.0

        #print('Cr-polar times 0:9 = ', timesecs_cr[0:9])
        #print('Cr-polar backscatter 0:4, 0:4 = ', vals_cr[0:4,0:4,0])
        #print('Cr-polar wind 0:4, 0:4 = ', vals_cr[0:4,0:4,1])
        #print('Cr-polar SNR+1 0:4, 0:4 = ', vals_cr[0:4,0:4,2])

        nn_start = 1    #Cr-polar
        new_ind = halo_check_monotonic(timesecs_cr, vals_cr, elev_angle_cr, n_cr, nn_start, datestring, ngates, errorfile, midnight_file_root)
        #print('cr-polar new_ind = ',new_ind)
        if new_ind >= 0:        #Non-monotonic data found
            n_cr = n_cr - 1
            timesecs_cr = timesecs_cr[0:new_ind]
            #azi_arr = azi_arr[0:new_ind,:]
            elev_angle_cr = elev_angle_cr[0:new_ind]
            vals_cr = vals_cr[0:new_ind,:,:]


    if n_co > 0:

        # ---------------------
        # Setup height variable
        # ---------------------
        h    = np.zeros((ngates))
        ind_range = np.arange(ngates)
        h[:] = range(ngates)
        h    = h+0.5	#To refer height to centre of bin, not bottom (or top)
        #h   *= 0.001*rgate*sin(elev_angle_co[0]*pi/180)	#Elev angle is saved every scan - should check it doesn't change through day
        h   *= rgate	#Don't scale by elev_angle, it's range, not height
        #print('Range values 0:9 = ', h[0:9])

        # ----------
        # Processing
        # ----------
        # Exclude time averaging code from this version. Can be added later
        tinc = 0	#Forces no time averaging
        if tinc < 10: 
            # -------------------
            # Omit time averaging
            # -------------------
            #print('Skipping time averaging')

            #Calculate valid min and max values for variable attributes
            bs_co_real = vals_co[:,:,0][np.where(vals_co[:,:,0] >= 0)]	#bs_co_real becomes 1D at this point, ok if just need values  
            bs_cr_real = vals_cr[:,:,0][np.where(vals_cr[:,:,0] >= 0)]  
            bs_vad_real = vals_vad[:,:,0][np.where(vals_vad[:,:,0] >= 0)]  
            #Setup qc_flag values, don't have any restrictions for SNR, so qc_flag_co, _cr and _vad are for backscattering coefficient
            qc_flag_co=np.ones((n_co,ngates))
            qc_flag_cr=np.ones((n_cr,ngates))
            qc_flag_vad=np.ones((n_vad,ngates))
            qc_flag_vel_co=np.ones((n_co,ngates))
            qc_flag_vel_cr=np.ones((n_cr,ngates))
            qc_flag_vel_vad=np.ones((n_vad,ngates))
            for kk in range(n_co):
                for ll in range(ngates):
                    if ll < (ngates-1):
                        if np.fabs((vals_co[kk,ll+1,1] - vals_co[kk,ll,1])) > 5.0:	#Unrealistic wind shear defined as 5 ms-1 change between 30m bins
                            qc_flag_vel_co[kk,ll] = 3 
                    if vals_co[kk,ll,1] > 19.0 or vals_co[kk,ll,1] < -19.0:
                        qc_flag_vel_co[kk,ll] = 2 
                    if vals_co[kk,ll,0] < 0:
                        qc_flag_co[kk,ll] = 3
                    if vals_co[kk,ll,0] < 1.0E-07 or vals_co[kk,ll,0] > 1000:
                        qc_flag_co[kk,ll] = 2
            beta_co_good = vals_co[:,:,0][np.where(qc_flag_co == 1)]
            vel_co_good = vals_co[:,:,1][np.where(qc_flag_vel_co == 1)]
            for kk in range(n_cr):
                for ll in range(ngates):
                    if ll < (ngates-1):
                        if np.fabs((vals_cr[kk,ll+1,1] - vals_cr[kk,ll,1])) > 5.0:	#Unrealistic wind shear defined as 5 ms-1 change between 30m bins
                            qc_flag_vel_cr[kk,ll] = 3 
                    if vals_cr[kk,ll,1] > 19.0 or vals_cr[kk,ll,1] < -19.0:
                        qc_flag_vel_cr[kk,ll] = 2 
                    if vals_cr[kk,ll,0] < 0:
                        qc_flag_cr[kk,ll] = 3
                    if vals_cr[kk,ll,0] < 1.0E-07 or vals_cr[kk,ll,0] > 1000:
                        qc_flag_cr[kk,ll] = 2
            beta_cr_good = vals_cr[:,:,0][np.where(qc_flag_cr == 1)]
            vel_cr_good = vals_cr[:,:,1][np.where(qc_flag_vel_cr == 1)]
            for kk in range(n_vad):
                for ll in range(ngates):
                    if ll < (ngates-1):
                        if np.fabs((vals_vad[kk,ll+1,1] - vals_vad[kk,ll,1])) > 5.0:	#Unrealistic wind shear defined as 5 ms-1 change between 30m bins
                            qc_flag_vel_vad[kk,ll] = 3 
                    if vals_vad[kk,ll,1] > 19.0 or vals_vad[kk,ll,1] < -19.0:
                        qc_flag_vel_vad[kk,ll] = 2 
                    if vals_vad[kk,ll,0] < 0:
                        qc_flag_vad[kk,ll] = 3
                    if vals_vad[kk,ll,0] < 1.0E-07 or vals_vad[kk,ll,0] > 1000:
                        qc_flag_vad[kk,ll] = 2
            beta_vad_good = vals_vad[:,:,0][np.where(qc_flag_vad == 1)]
            vel_vad_good = vals_vad[:,:,1][np.where(qc_flag_vel_vad == 1)]

            #print('min, max co = ', np.amin(beta_co_good), np.amax(beta_co_good), np.amin(vel_co_good), np.amax(vel_co_good))
            #print('min, max cr = ', np.amin(beta_cr_good), np.amax(beta_cr_good), np.amin(vel_cr_good), np.amax(vel_cr_good))
            #print('min, max VAD = ', np.amin(beta_vad_good), np.amax(beta_vad_good), np.amin(vel_vad_good), np.amax(vel_vad_good))


            sampling_interval = str(timesecs_co[1]-timesecs_co[0])+' seconds'

        if infile_type == 2:	#If writing a netCDF file with cleaned data
            vals_co_cleaned = missing_value * np.ones((n_co,ngates,3))  #Default data array, will be unchanged if file doesn't exist
            # ----------------------------------
            # Open cleaned netCDF file to read 
            # ----------------------------------
            cleaned_nc_file = cleaned_nc_path_in +  str(year) + '/' + datestring + '_chilbolton_halo-doppler-lidar-118-stare-co.nc'
            if os.path.isfile(cleaned_nc_file):
                print('Opening cleaned NetCDF file ' + cleaned_nc_file)
                ncfile=nc4.Dataset(cleaned_nc_file, 'r', format='NETCDF3_CLASSIC')
                temp_time = ncfile.variables['time']
                temp_beta = ncfile.variables['beta']
                temp_vel = ncfile.variables['v']        #Think this is the correct one to read, with snr=1 filtering?
                temp_sig = ncfile.variables['signal']
                temp_beta_error = ncfile.variables['beta_error']
                temp_vel_error = ncfile.variables['v_error']
                cleaned_time = temp_time[:]     #Otherwise it disappears when you close the .nc file
                cleaned_beta = cal_factor * temp_beta[:]        #Calibrate cleaned data here
                cleaned_vel = temp_vel[:]
                cleaned_snr = temp_sig[:]
                fr_beta_error = temp_beta_error[:]      #Defined in file as fractional error (no units)
                fr_vel_error = temp_vel_error[:]        #Defined in file as velocity error (ms-1)
                cleaned_time *= 3600.0
                fr_vel_error = fr_vel_error/cleaned_vel #Convert velocity error to fractional
                print('Closing cleaned NetCDF file')
                ncfile.close()

                n_co_cleaned_in = len(cleaned_time)     #Save to check it matches uncleaned data
                print('Number of raw, cleaned co-polar profiles = ', n_co, n_co_cleaned_in)
                #Check for differences in the number of time values between this processing and cleaned files

                if n_co_cleaned_in == n_co:     #Easy case - am assuming that the timestamps of cleaned file match this processing
                    vals_co_cleaned[:,:,0] = cleaned_beta
                    vals_co_cleaned[:,:,1] = cleaned_vel
                    vals_co_cleaned[:,:,2] = cleaned_snr
                    co_cleaned_fail = 0
                    #print('vals_co_cleaned at 1640, 30-40 = ', vals_co_cleaned[1640, 30:40, 0])

                if (n_co_cleaned_in - n_co) == 1:       #More points in cleaned file, e.g. 20201206, where cleaned file erroneously has point 1s into 20201207 wrongly assigned to the start of 20201207
                    matching_index = np.argmin(np.absolute(cleaned_time - timesecs_co[0]))      #Are times in same units?
                    #print('n_co_cleaned_in - n_co = 1 matching_index = ', matching_index)

                    vals_co_cleaned[:,:,0] = cleaned_beta[matching_index:n_co_cleaned_in,:]
                    vals_co_cleaned[:,:,1] = cleaned_vel[matching_index:n_co_cleaned_in,:]
                    vals_co_cleaned[:,:,2] = cleaned_snr[matching_index:n_co_cleaned_in,:]
                    fr_beta_error = fr_beta_error[matching_index:n_co_cleaned_in,:]
                    fr_vel_error = fr_vel_error[matching_index:n_co_cleaned_in,:]

                    err_file=open(errorfile,'a')        #Error log text file
                    err_file.write(oformaterr1 % (datestring))
                    err_file.write("%s" % "  0 ")
                    err_file.write(oformaterr3 % (n_co))
                    err_file.write("%s" % " ")
                    err_file.write(oformaterr3 % (n_co_cleaned_in))
                    err_file.write("%s" % " ")
                    err_file.write(oformatmid1 % (cleaned_time[0]))     #Seconds
                    err_file.write("%s" % " Co-cleaned\n")
                    err_file.close()

                    co_cleaned_fail = 0

                if (n_co_cleaned_in - n_co) == -1:      #Less points in cleaned file, e.g. 20201207, where cleaned file erroneously has point missing 1s after start
                    matching_index = np.argmin(np.absolute(timesecs_co - cleaned_time[0]))      #Are times in same units?
                    #print('n_co_cleaned_in - n_co = -1 matching_index = ', matching_index)

                    vals_co_cleaned[matching_index:(n_co + matching_index - 1),:,0] = cleaned_beta
                    vals_co_cleaned[matching_index:(n_co + matching_index - 1),:,1] = cleaned_vel
                    vals_co_cleaned[matching_index:(n_co + matching_index - 1),:,2] = cleaned_snr
                    #fr_beta_error[0,:] = 1.0    #Set first time of error array to 1.0
                    #fr_vel_error[0,:] = 1.0     #Set first time of error array to 1.0
                    fr_error_temp = np.ones((n_co, ngates))
                    fr_error_temp[matching_index:(n_co + matching_index - 1),:] =fr_beta_error
                    fr_beta_error = fr_error_temp
                    fr_error_temp = np.ones((n_co, ngates))
                    fr_error_temp[matching_index:(n_co + matching_index - 1),:] =fr_vel_error
                    fr_vel_error = fr_error_temp

                    err_file=open(errorfile,'a')        #Error log text file
                    err_file.write(oformaterr1 % (datestring))
                    err_file.write("%s" % "  0 ")
                    err_file.write(oformaterr3 % (n_co))
                    err_file.write("%s" % " ")
                    err_file.write(oformaterr3 % (n_co_cleaned_in))
                    err_file.write("%s" % " ")
                    err_file.write(oformatmid1 % (timesecs_co[0]))      #Seconds
                    err_file.write("%s" % " Co-cleaned\n")
                    err_file.close()

                    co_cleaned_fail = 1

                #If there is more than 1 data point difference between uncleaned and clean, make all cleaned data missing
                if np.absolute(n_co_cleaned_in - n_co) > 1:
                    co_cleaned_fail = 2

            else:       #Cleaned file not found
                co_cleaned_fail = 2

            #Setup cleaned_qc_flag values
            #Do we want to qc_flag = keep 3 as below inst threshold? I think it's no longer relevant
            qc_flag_co_cleaned=np.ones((n_co,ngates))   #Have assumed no points from cleaned file is same as we have read for raw data
            qc_flag_vel_co_cleaned=np.ones((n_co,ngates))       #Have assumed no points from cleaned file is same as we have read for raw data
            for kk in range(n_co):
                for ll in range(ngates):
                    if ll < (ngates-1):
                        if np.fabs((vals_co_cleaned[kk,ll+1,1] - vals_co_cleaned[kk,ll,1])) > 5.0:      #Unrealistic wind shear defined as 5 ms-1 change between 30m bins
                            qc_flag_vel_co_cleaned[kk,ll] = 3
                    if vals_co_cleaned[kk,ll,1] > 19.0 or vals_co_cleaned[kk,ll,1] < -19.0:
                        qc_flag_vel_co_cleaned[kk,ll] = 2
                    #if vals_co[kk,ll,0] < 0:
                        #qc_flag_co[kk,ll] = 3
                    if vals_co_cleaned[kk,ll,0] < 1.0E-07 or vals_co_cleaned[kk,ll,0] > 1000:
                        qc_flag_co_cleaned[kk,ll] = 2
            #Add quality flag based on fractional error in data, as read from cleaned file
            #By putting it here, it will over-ride the previous flags of 2 (and 3 for velocity data)
            #I feel this is fair, because some of those flags will be caused by noisy data
            qc_flag_co_cleaned = np.where(np.absolute(fr_beta_error) <= uncertainty_threshold, qc_flag_co_cleaned, 3)
            qc_flag_vel_co_cleaned = np.where(np.absolute(fr_vel_error) <= uncertainty_threshold, qc_flag_vel_co_cleaned, 4)
            n_speckles_co = 0
            n_speckles_co_vel = 0
            #Look for speckles - points with fr_error less than threshold, but surrounded by at least 3 out of 4 points with fr_error greater than threshold
            for kk in range(n_co):
                for ll in range(ngates):
                    if kk > 0 and kk < (n_co-1) and ll > 0 and ll < (ngates-1):
                        if qc_flag_co_cleaned[kk,ll] == 1:	#Not flagged as bad already
                            array_to_check = np.asarray([qc_flag_co_cleaned[kk-1,ll], qc_flag_co_cleaned[kk+1,ll], qc_flag_co_cleaned[kk,ll-1], qc_flag_co_cleaned[kk,ll+1]])
                            yes_no_check =np.where(array_to_check == 3, 1, 0)	#Check whether each of the surrounding points has a fractional error above the threshold
                            if round(np.sum(yes_no_check)) >= 3:    #If at least 3 points have a error above the threshold
                                qc_flag_co_cleaned[kk,ll] = 4
                                n_speckles_co = n_speckles_co + 1 
                        if qc_flag_vel_co_cleaned[kk,ll] == 1:	#Not flagged as bad already
                            array_to_check = np.asarray([qc_flag_vel_co_cleaned[kk-1,ll], qc_flag_vel_co_cleaned[kk+1,ll], qc_flag_vel_co_cleaned[kk,ll-1], qc_flag_vel_co_cleaned[kk,ll+1]])
                            yes_no_check =np.where(array_to_check == 4, 1, 0)	#Check whether each of the surrounding points has a fractional error above the threshold
                            if round(np.sum(yes_no_check)) >= 3:    #If at least 3 points have a error above the threshold
                                qc_flag_vel_co_cleaned[kk,ll] = 5 
                                n_speckles_co_vel = n_speckles_co_vel + 1 
            print('Number of backscatter, velocity co-polar points identified as speckles = ', n_speckles_co, n_speckles_co_vel)
            #Change these to 5,6 (with speckle flag) from 4,5 (without speckle flag) 
            if co_cleaned_fail == 1:
                qc_flag_co_cleaned[0,:] = 5
                qc_flag_vel_co_cleaned[0,:] = 6
            if co_cleaned_fail == 2:
                qc_flag_co_cleaned[:,:] = 5
                qc_flag_vel_co_cleaned[:,:] = 6
            #Subset of qc_flag = 1 values to determine valid min, max
            beta_co_cleaned_good = vals_co_cleaned[:,:,0][np.where(qc_flag_co_cleaned == 1)]
            vel_co_cleaned_good = vals_co_cleaned[:,:,1][np.where(qc_flag_vel_co_cleaned == 1)]
            #If data were all missing, insert 1 missing value to these cleaned arrays in order to be able to write file
            if len(beta_co_cleaned_good) == 0:
                beta_co_cleaned_good = missing_value
            if len(vel_co_cleaned_good) == 0:
                vel_co_cleaned_good = missing_value

            #--------------------------------------------
            #Read cleaned VAD files - note 1 per VAD scan
            #--------------------------------------------

            cleaned_vad_nc_path_in = "/data/lidar-doppler-halo/netcdf/calibrated/halo-doppler-lidar-118-vad75/" + str(year) + '/' + datestring + '/'
            files_cleaned_vad = os.listdir(cleaned_vad_nc_path_in)
            expr = re.compile('nc')
            infiles_cleaned_vad = [elem for elem in files_cleaned_vad if expr.search(elem)]
            infiles_cleaned_vad.sort()
            #print(infiles_cleaned_vad)
            nfiles_cleaned_vad = len(infiles_cleaned_vad)
            if nfiles_cleaned_vad == 0:
                print('No cleaned VAD files to process')

            n_vad_cleaned = 0

            if nfiles_cleaned_vad > 0:  #Some cleaned VAD files were found

                time_vad_cleaned = np.zeros(n_vad)
                vals_vad_cleaned_in = np.zeros((n_vad,ngates,3))
                fr_beta_error_vad_in = np.zeros((n_vad,ngates))
                fr_vel_error_vad_in = np.zeros((n_vad,ngates))
                for nf in range(nfiles_cleaned_vad):
                    ncfile=nc4.Dataset((cleaned_vad_nc_path_in + infiles_cleaned_vad[nf]), 'r', format='NETCDF3_CLASSIC')
                    temp_time = ncfile.variables['time']
                    temp_beta = ncfile.variables['beta']
                    temp_vel = ncfile.variables['v']    #Think this is the correct one to read, with snr=1 filtering?
                    temp_sig = ncfile.variables['signal']
                    temp_beta_error = ncfile.variables['beta_error']
                    temp_vel_error = ncfile.variables['v_error']
                    cleaned_time = temp_time[:] #Not used with VAD data, so not converted from hours to seconds
                    cleaned_beta = cal_factor * temp_beta[:]    #Apply calibration here
                    cleaned_vel = temp_vel[:]
                    cleaned_snr = temp_sig[:]
                    fr_beta_error = temp_beta_error[:]      #Defined in file as fractional error (no units)
                    fr_vel_error = temp_vel_error[:]        #Defined in file as velocity error (ms-1)
                    #print 'clean_beta shape = ',cleaned_beta.shape
                    ncfile.close()
                    n_this_vad = len(cleaned_time)
                    time_vad_cleaned[n_vad_cleaned:(n_vad_cleaned + n_this_vad)] = cleaned_time
                    vals_vad_cleaned_in[n_vad_cleaned:(n_vad_cleaned + n_this_vad),:,0] = cleaned_beta
                    vals_vad_cleaned_in[n_vad_cleaned:(n_vad_cleaned + n_this_vad),:,1] = cleaned_vel
                    vals_vad_cleaned_in[n_vad_cleaned:(n_vad_cleaned + n_this_vad),:,2] = cleaned_snr
                    fr_beta_error_vad_in[n_vad_cleaned:(n_vad_cleaned + n_this_vad),:] = fr_beta_error 
                    fr_vel_error_vad_in[n_vad_cleaned:(n_vad_cleaned + n_this_vad),:] = fr_beta_error 

                    n_vad_cleaned = n_vad_cleaned + n_this_vad

                print('Number of raw, cleaned VAD profiles = ',n_vad, n_vad_cleaned)

                #Check if cleaned VAD time is monotonic - if will have issues writing data to the correct uncleaned times from .proc files
                time_vad_diff = np.diff(time_vad_cleaned)
                if np.amin(time_vad_diff) < 0:  #So there's at least 1 point with a negative difference
                    no_time_discrep = np.sum(np.where(time_vad_diff < 0, 1, 0))
                    #print('Cleaned VAD time not monotonic, number of time discrepancies = ', no_time_discrep)
                    if no_time_discrep == 1:
                        vals_vad_cleaned = missing_value * np.ones((n_vad,ngates,3))
                        fr_beta_error_vad = missing_value * np.ones((n_vad,ngates))
                        fr_vel_error_vad = missing_value * np.ones((n_vad,ngates))
                        vad_diff_neg = np.nonzero(time_vad_diff < 0)
                        ind_vad_diff_neg = vad_diff_neg[0][0] + 1
                        #print('Cleaned VAD time not monotonic, offset at index, time range = ', ind_vad_diff_neg, time_vad_cleaned[(ind_vad_diff_neg-1):(ind_vad_diff_neg+2)])
                        if ind_vad_diff_neg <= 24:      #Discrepancy within 24 points of start of file (i.e. 1 VAD scan)
                            #Take out data from the input arrays up to the point where time becomes monotonic
                            vals_vad_cleaned_in = vals_vad_cleaned_in[ind_vad_diff_neg:,:,:]
                            fr_beta_error_vad_in = fr_beta_error_vad_in[ind_vad_diff_neg:,:]
                            fr_vel_error_vad_in = fr_vel_error_vad_in[ind_vad_diff_neg:,:]
                            vals_vad_cleaned[0:(n_vad-ind_vad_diff_neg),:,:] = vals_vad_cleaned_in
                            fr_beta_error_vad[0:(n_vad-ind_vad_diff_neg),:] = fr_beta_error_vad_in 
                            fr_vel_error_vad[0:(n_vad-ind_vad_diff_neg),:] = fr_vel_error_vad_in 
                            non_mono_qc = 1

                            err_file=open(errorfile,'a')        #Error log text file
                            err_file.write(oformaterr1 % (datestring))
                            err_file.write("%s" % "  2 ")
                            err_file.write(oformaterr3 % (n_vad))
                            err_file.write("%s" % " ")
                            err_file.write(oformaterr3 % (n_vad_cleaned))
                            err_file.write("%s" % " ")
                            err_file.write(oformatmid1 % (3600.0*time_vad_cleaned[0]))  #Seconds
                            err_file.write("%s" % " VAD-cleaned\n")
                            err_file.close()
                        else:   #After 24 points, which means it's hard to be sure how profiles are ordered through day
                            print('Non-monotonic time index > 24 (1 scan), will not try to match cleaned VAD data to times')
                            #Set an indicator that sets the qc_flag to a value showing this
                            non_mono_qc = 2

                    else:
                        print('More than 1 monotonic point detected, will not try to match cleaned VAD data to times')
                        vals_vad_cleaned = missing_value * np.ones((n_vad,ngates,3))
                        #Set an indicator that sets the qc_flag to a value showing this
                        non_mono_qc = 2
                else:   #Data are monotonic
                    #Check that uncleaned and cleaned VAD files have same number of points
                    if n_vad == n_vad_cleaned:
                        vals_vad_cleaned = vals_vad_cleaned_in
                        fr_beta_error_vad = fr_beta_error_vad_in 
                        fr_vel_error_vad = fr_vel_error_vad_in 
                        non_mono_qc = 0
                    else:       #
                        print('Unequal numbers of uncleaned and cleaned VAD scans, will not try to match cleaned VAD data to times')
                        vals_vad_cleaned = missing_value * np.ones((n_vad,ngates,3))
                        #Set an indicator that sets the qc_flag to a value showing this
                        non_mono_qc = 2

                #Setup cleaned_qc_flag values for VAD
                #Do we want to qc_flag = keep 3 as below inst threshold? I think it's no longer relevant
                qc_flag_vad_cleaned=np.ones((n_vad,ngates))
                qc_flag_vel_vad_cleaned=np.ones((n_vad,ngates))
                for kk in range(n_vad):
                    for ll in range(ngates):
                        if ll < (ngates-1):
                            if np.fabs((vals_vad_cleaned[kk,ll+1,1] - vals_vad_cleaned[kk,ll,1])) > 5.0:    #Unrealistic wind shear defined as 5 ms-1 change between 30m bins
                                qc_flag_vel_vad_cleaned[kk,ll] = 3
                        if vals_vad_cleaned[kk,ll,1] > 19.0 or vals_vad_cleaned[kk,ll,1] < -19.0:
                            qc_flag_vel_vad_cleaned[kk,ll] = 2
                        #if vals_co[kk,ll,0] < 0:
                            #qc_flag_co[kk,ll] = 3
                        if vals_vad_cleaned[kk,ll,0] < 1.0E-07 or vals_vad_cleaned[kk,ll,0] > 1000:
                            qc_flag_vad_cleaned[kk,ll] = 2

                #Add quality flag based on fractional error in data, as read from cleaned file
                #By putting it here, it will over-ride the previous flags of 2 (and 3 for velocity data)
                #I feel this is fair, because some of those flags will be caused by noisy data
                qc_flag_vad_cleaned = np.where(np.absolute(fr_beta_error_vad) <= uncertainty_threshold, qc_flag_vad_cleaned, 3)
                qc_flag_vel_vad_cleaned = np.where(np.absolute(fr_vel_error_vad) <= uncertainty_threshold, qc_flag_vel_vad_cleaned, 4)

                #Look for speckles - points with fr_error less than threshold, but surrounded by at least 3 out of 4 points with fr_error greater than threshold
                n_speckles_vad = 0
                n_speckles_vad_vel = 0
                for kk in range(n_vad):
                    for ll in range(ngates):
                        if kk > 0 and kk < (n_vad-1) and ll > 0 and ll < (ngates-1):
                            if qc_flag_vad_cleaned[kk,ll] == 1:      #Not flagged as bad already
                                array_to_check = np.asarray([qc_flag_vad_cleaned[kk-1,ll], qc_flag_vad_cleaned[kk+1,ll], qc_flag_vad_cleaned[kk,ll-1], qc_flag_vad_cleaned[kk,ll+1]])
                                yes_no_check =np.where(array_to_check == 3, 1, 0)   #Check whether each of the surrounding points has a fractional error above the threshold
                                if round(np.sum(yes_no_check)) >= 3:    #If at least 3 points have a error above the threshold
                                    qc_flag_vad_cleaned[kk,ll] = 4
                                    n_speckles_vad = n_speckles_vad + 1
                            if qc_flag_vel_vad_cleaned[kk,ll] == 1:  #Not flagged as bad already
                                array_to_check = np.asarray([qc_flag_vel_vad_cleaned[kk-1,ll], qc_flag_vel_vad_cleaned[kk+1,ll], qc_flag_vel_vad_cleaned[kk,ll-1], qc_flag_vel_vad_cleaned[kk,ll+1]])
                                yes_no_check =np.where(array_to_check == 4, 1, 0)   #Check whether each of the surrounding points has a fractional error above the threshold
                                if round(np.sum(yes_no_check)) >= 3:    #If at least 3 points have a error above the threshold
                                    qc_flag_vel_vad_cleaned[kk,ll] = 5
                                    n_speckles_vad_vel = n_speckles_vad_vel + 1
                print('Number of backscatter, velocity VAD co-polar points identified as speckles = ', n_speckles_vad, n_speckles_vad_vel)

                if non_mono_qc == 1:
                    qc_flag_vad_cleaned[(n_vad-ind_vad_diff_neg):,:] = 5
                    qc_flag_vel_vad_cleaned[(n_vad-ind_vad_diff_neg):,:] = 6
                if non_mono_qc == 2:
                    qc_flag_vad_cleaned[:,:] = 5
                    qc_flag_vel_vad_cleaned[:,:] = 6

                #Subset of qc_flag = 1 values to determine valid min, max
                beta_vad_cleaned_good = vals_vad_cleaned[:,:,0][np.where(qc_flag_vad_cleaned == 1)]
                vel_vad_cleaned_good = vals_vad_cleaned[:,:,1][np.where(qc_flag_vel_vad_cleaned == 1)]
                #If data were all missing, insert 1 missing value to these cleaned arrays in order to be able to write file
                if len(beta_vad_cleaned_good) == 0:
                    beta_vad_cleaned_good = missing_value
                if len(vel_vad_cleaned_good) == 0:
                    vel_vad_cleaned_good = missing_value


            else:       #No cleaned VAD files found
                vals_vad_cleaned = missing_value * np.ones((n_vad,ngates,3))
                #Set an indicator that sets the qc_flag to a value showing this
                non_mono_qc = 2


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        for nn in range(3):	# Co-polar, cross-polar, VAD
        #for nn in range(2):	# Co-polar, cross-polar. VAD taken out as files are big
            if nn == 0:
                if n_co > 0:
                    ok_to_write[0,nn] = 1
                no_t_av = n_co
                time_av = timesecs_co
                elev_av = elev_angle_co
                elev_arr=np.transpose(np.tile(elev_av,(1,1)))	#For AMF requirements, h has time dimension, same shape as bscatt values
                azi_arr = np.zeros((no_t_av,1))
                vals_av = vals_co
                qc_flag = qc_flag_co
                qc_flag_vel = qc_flag_vel_co
                beta_good = beta_co_good
                vel_good = vel_co_good
                no_files_to_process = infile_type
                if infile_type == 2 > 0:	#Cleaned data
                    if n_co_cleaned_in > 0:
                        ok_to_write[1,nn] = 1
                        vals_cleaned = vals_co_cleaned
                        #print('vals_cleaned at 1640, 30-40 = ', vals_cleaned[1640, 30:40, 0])
                        beta_cleaned_good = beta_co_cleaned_good
                        vel_cleaned_good = vel_co_cleaned_good
                        snr_cleaned_good = vals_co_cleaned[:,:,2]
                        qc_flag_cleaned = qc_flag_co_cleaned
                        qc_flag_vel_cleaned = qc_flag_vel_co_cleaned
                data_name = ["_aerosol-backscatter-radial-winds_fixed_co_standard_v1.0","_aerosol-backscatter-radial-winds_fixed_co_advanced_v1.0"]
                template_file_path_lv1 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-co-lv1_template.yaml"	#YAML template
                template_file_path_lv2 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-co-lv2_template.yaml"	#YAML template

            elif nn == 1 and n_cr > 0:
                if n_cr > 0:
                    ok_to_write[0,nn] = 1
                no_t_av = n_cr
                time_av = timesecs_cr
                elev_av = elev_angle_cr
                elev_arr=np.transpose(np.tile(elev_av,(1,1)))	#For AMF requirements, h has time dimension, same shape as bscatt values
                azi_arr = np.zeros((no_t_av,1))
                vals_av = vals_cr
                qc_flag = qc_flag_cr
                qc_flag_vel = qc_flag_vel_cr
                beta_good = beta_cr_good
                vel_good = vel_cr_good
                no_files_to_process = 1	#Only produce standard data from cross-polar, there is no cleaned
                data_name = ["_aerosol-backscatter-radial-winds_fixed_cr_standard_v1.0","_aerosol-backscatter-radial-winds_fixed_cr_advanced_v1.0"]	#As yet, no advanced cross polar data, so 2nd name is ignored
                template_file_path_lv1 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-cr-lv1_template.yaml"	#YAML template
                template_file_path_lv2 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-cr-lv1_template.yaml"	#YAML template, not used at present
                #template_file_path_lv2 = "/home/jla/python/halo/cfarr-lidar-halo-lv2_metadata-template.yaml"	#YAML template, not used at present

            else: #nn == 2:
                if n_vad > 0:
                    ok_to_write[0,nn] = 1
                nangles = 24
                no_t_av = n_vad
                time_av = timesecs_vad
                elev_arr = elev_angle_vad
                azi_arr = azi_angle_vad
                vals_av = vals_vad
                qc_flag = qc_flag_vad
                qc_flag_vel = qc_flag_vel_vad
                beta_good = beta_vad_good
                vel_good = vel_vad_good
                no_files_to_process = infile_type 
                if infile_type == 2:	#Cleaned data
                    if n_vad_cleaned > 0:
                        ok_to_write[1,nn] = 1
                        vals_cleaned = vals_vad_cleaned
                        beta_cleaned_good = beta_vad_cleaned_good
                        vel_cleaned_good = vel_vad_cleaned_good
                        snr_cleaned_good = vals_vad_cleaned[:,:,2]
                        qc_flag_cleaned = qc_flag_co_cleaned
                        qc_flag_cleaned = qc_flag_vad_cleaned
                        qc_flag_vel_cleaned = qc_flag_vel_vad_cleaned
                data_name = ["_aerosol-backscatter-radial-winds_vad_standard_v1.0","_aerosol-backscatter-radial-winds_vad_advanced_v1.0"]	#As yet, no advanced cross polar data, so 2nd name is ignored
                template_file_path_lv1 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-co-lv1-vad_template.yaml"	#YAML template
                template_file_path_lv2 = "/home/chilbolton_software/python/ncas_python/halo/ncas-lidar-dop-3-co-lv2-vad_template.yaml"	#YAML template
                #template_file_path_lv1 = "/home/jla/python/halo/cfarr-lidar-halo-co-lv1-vad_metadata-template.yaml"	#YAML template
                #template_file_path_lv2 = "/home/jla/python/halo/cfarr-lidar-halo-co-lv2-vad_metadata-template.yaml"	#YAML template
	

            tt = (year,month,day,0,0,0,0,0,0)
            tt_upper = (year,month,(day+1),0,0,0,0,0,0)
            day_start_epoch = time.mktime(tt)
            #t_valid_min = day_start_epoch
            #t_valid_max = time.mktime(tt_upper) 
            yr_arr=np.zeros(no_t_av)+year
            mn_arr=np.zeros(no_t_av)+month
            dy_arr=np.zeros(no_t_av)+day
            dyyr_arr=np.zeros(no_t_av)
            epoch_timesecs=np.zeros(no_t_av)
            hr_arr=np.zeros(no_t_av,dtype=np.int32)
            mi_arr=np.zeros(no_t_av,dtype=np.int32)
            sc_arr=np.zeros(no_t_av,dtype=np.float32)
            epoch_timesecs = day_start_epoch + time_av
            for mm in range(no_t_av):
                #epoch_timesecs[mm] = day_start_epoch + time_av[mm]
                tt = time.gmtime(epoch_timesecs[mm])
                hr_arr[mm] = np.int32(tt[3])
                mi_arr[mm] = np.int32(tt[4])
                sc_arr[mm] = time_av[mm]-3600*hr_arr[mm]-60*mi_arr[mm]	#Keeps digits after decimal place, unlike tt[5] 
                dyyr_arr[mm] = np.float32(tt[7]+time_av[mm]/86400.0)
                if time_av[mm] >= 86400.0:
                    sc_arr[mm] = sc_arr[mm] - 86400.0	#Next day, should only be present for VADs
                    dy_arr[mm] = dy_arr[mm] + 1	#Value only based on "day" calculated from start date
                    dyyr_arr[mm] = dyyr_arr[mm] - 1.0	#Value based on tt[7] which accounts for new day 


            #h_arr=np.tile(h,(time_av.size,1))	#For AMF requirements, h has time dimension, same shape as bscatt values
            ##h_arr=np.transpose(h_arr)
            #h_arr=broadcast(h,no_t_av,ngates)

            h_arr = my_broadcast(h,no_t_av,ngates,nangles)
            #h_arr_int1=np.broadcast_to(h,(nangles,time_av.size,ngates))
            #h_arr_int2=np.swapaxes(h_arr_int1,0,2)
            #h_arr=np.swapaxes(h_arr_int2,0,1)

            #------------Setting up YAML data object----------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules

            tfp = [template_file_path_lv1, template_file_path_lv2]

            for n_outfile in range(no_files_to_process):	#Standard and cleaned loop

                if ok_to_write[n_outfile, nn] > 0: 
                    handler = module_data_object_python3.Handler()
                    data_object_id = handler.load_template(tfp[n_outfile])
                    print('Data object ID = ', data_object_id)
                    #handler.show_requirements_for_template(data_object_id)

                    # -------------------------------
                    # Open new netCDF file for output
                    # -------------------------------
                    out_file   = os.path.join(path_out,cfarr_head+datestring+data_name[n_outfile]+'.nc')
                    print('Opening new NetCDF file ' + out_file)

                    lengths_of_dimensions = {"time": no_t_av, "index_of_range": ngates, "index_of_angle": nangles}
                    uname              = os.uname()
                    nodename           = uname[1]
                    os_name            = uname[0]
                    os_release         = uname[2]
                    computer_id        = nodename + ' under ' + os_name + ' ' + os_release
                    user_id            = pwd.getpwuid(os.getuid())[0]


                    if n_outfile == 0:	#Not cleaned	    
                        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],int(sc_arr[0])), "end_time": datetime.datetime(year,month,day,hr_arr[no_t_av-1],mi_arr[no_t_av-1],int(sc_arr[no_t_av-1])), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "h_min": np.amin(h_arr), "h_max": np.amax(h_arr), "zen_min": np.amin(elev_arr), "zen_max": np.amax(elev_arr), "azi_min": np.amin(azi_arr), "azi_max": np.amax(azi_arr), "min_bs": np.amin(beta_good), "max_bs": np.amax(beta_good), "min_v": np.amin(vel_good), "max_v": np.amax(vel_good), "min_snr": np.amin(vals_av[:,:,2]), "max_snr": np.amax(vals_av[:,:,2]), "uname": user_id, "codename": __file__, "machinename": computer_id}

                    else:	#Cleaned
                        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],int(sc_arr[0])), "end_time": datetime.datetime(year,month,day,hr_arr[no_t_av-1],mi_arr[no_t_av-1],int(sc_arr[no_t_av-1])), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "h_min": np.amin(h_arr), "h_max": np.amax(h_arr), "zen_min": np.amin(elev_arr), "zen_max": np.amax(elev_arr), "azi_min": np.amin(azi_arr), "azi_max": np.amax(azi_arr), "min_clean_bs": np.amin(beta_cleaned_good), "max_clean_bs": np.amax(beta_cleaned_good), "min_clean_vel": np.amin(vel_cleaned_good), "max_clean_vel": np.amax(vel_cleaned_good), "min_clean_snr": np.amin(snr_cleaned_good), "max_clean_snr": np.amax(snr_cleaned_good), "uname": user_id, "codename": __file__, "machinename": computer_id, "bs_cal": str(cal_factor)}
 
                    data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
                    data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
                    data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]
                    data_object["variables"]["year"]["values"][:] = yr_arr[:]
                    data_object["variables"]["month"]["values"][:] = mn_arr[:]
                    data_object["variables"]["day"]["values"][:] = dy_arr[:]
                    data_object["variables"]["hour"]["values"][:] = hr_arr[:]
                    data_object["variables"]["minute"]["values"][:] = mi_arr[:]
                    data_object["variables"]["second"]["values"][:] = sc_arr[:]
                    #data_object["variables"]["range"]["values"][:] = h[:]	#When range is 1D
                    data_object["variables"]["range"]["values"][:] = h_arr[:]	#When range is 2D
                    data_object["variables"]["index_of_range"]["values"][:] = ind_range 
                    data_object["variables"]["sensor_view_angle_instrument_frame"]["values"][:] = elev_arr[:] 
                    data_object["variables"]["sensor_view_angle_earth_frame"]["values"][:] = elev_arr[:] 
                    data_object["variables"]["sensor_azimuth_angle_instrument_frame"]["values"][:] = azi_arr[:]
                    data_object["variables"]["sensor_azimuth_angle_earth_frame"]["values"][:] = azi_arr[:]
                    #print('azimuth angle shape = ',np.zeros((no_t_av,1)).shape) 

                    if n_outfile == 0:	#Standard, no clean
                        data_arr_out = my_broadcast(vals_av[:,:,0],no_t_av,ngates,nangles)
                        data_object["variables"]["attenuated_aerosol_backscatter_coefficient"]["values"][:] = data_arr_out
                        data_arr_out = my_broadcast(vals_av[:,:,1],no_t_av,ngates,nangles)
                        data_object["variables"]["radial_velocity_of_scatterers_away_from_instrument"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(vals_av[:,:,2],no_t_av,ngates,nangles)
                        data_object["variables"]["signal_to_noise_ratio_plus_1"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(qc_flag,no_t_av,ngates,nangles)
                        data_object["variables"]["qc_flag_attenuated_aerosol_backscatter_coefficient"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(qc_flag_vel,no_t_av,ngates,nangles)
                        data_object["variables"]["qc_flag_radial_velocity_of_scatterers_away_from_instrument"]["values"][:] = data_arr_out 
                    else:	#FMI cleaned
                        data_arr_out = my_broadcast(vals_cleaned[:,:,0],no_t_av,ngates,nangles)
                        data_object["variables"]["attenuated_aerosol_backscatter_coefficient"]["values"][:] = data_arr_out
                        data_arr_out = my_broadcast(vals_cleaned[:,:,1],no_t_av,ngates,nangles)
                        data_object["variables"]["radial_velocity_of_scatterers_away_from_instrument"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(vals_cleaned[:,:,2],no_t_av,ngates,nangles)
                        data_object["variables"]["signal_to_noise_ratio_plus_1"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(qc_flag_cleaned,no_t_av,ngates,nangles)
                        data_object["variables"]["qc_flag_attenuated_aerosol_backscatter_coefficient"]["values"][:] = data_arr_out 
                        data_arr_out = my_broadcast(qc_flag_vel_cleaned,no_t_av,ngates,nangles)
                        data_object["variables"]["qc_flag_radial_velocity_of_scatterers_away_from_instrument"]["values"][:] = data_arr_out 

                    exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
                    print('Data object exit code = ', exit_code)

                    print('Changing permissions of file ',out_file)
                    #oscommand = "chmod g+w " + out_file
                    oscommand = "chmod 775 " + out_file
                    os.system(oscommand)


                    ##if nn == 0:     #Only plot co-polar files for now
                    #if n_outfile == 0:
                        #lidar_plot.lidar_plotting(nday, 2, 1)   #Plot raw data
                    #if n_outfile == 1:
                        #lidar_plot.lidar_plotting(nday, 2, 2)   #Plot cleaned data
                    #print(" ")

    print(" ")
    for n_outfile in range(no_files_to_process):
        lidar_plot.lidar_plotting(nday, 2, (n_outfile + 1))   #Plot raw data
        print(" ")

    nprof_file=open(nproffile,'a')	#Error log text file
    nprof_file.write(oformaterr1 % (datestring))
    nprof_file.write("%s" % " ")
    nprof_file.write(oformaterr3 % (n_co))
    nprof_file.write("%s" % " ")
    nprof_file.write(oformaterr3 % (n_cr))
    nprof_file.write("%s" % " ")
    nprof_file.write(oformaterr3 % (n_vad))
    if infile_type == 2:
        nprof_file.write("%s" % " ")
        nprof_file.write(oformaterr3 % (n_co_cleaned_in))
        nprof_file.write("%s" % " ")
        nprof_file.write(oformaterr3 % (n_vad_cleaned))
    nprof_file.write("%s" % "\n")
    nprof_file.close()


# ------------------------------------------------------------------------
# Define function to check whether Halo data are monotonic
# ------------------------------------------------------------------------
def halo_check_monotonic(time_mono, vals_mono, elev_mono, nvals_mono, nn_mono, datestring, ngates, errorfile, midnight_file_root):
#Produces a time index at which non-monotonic time was detected

    #---------------------------------------------------------------------
    #Check that time_av is monotonic
    #For fixed scans, save points after midnight to put in next day's file
    #For VADs, add 86400 seconds to time and keep whole VAD in this file
    #---------------------------------------------------------------------

    oformaterr1="%s"
    oformaterr2="%2i"
    oformaterr3="%5i"
    oformatmid1="%13.6e"

    time_diff = np.diff(time_mono)
    if np.amin(time_diff) < 0:  #So there's a point with a negative difference
        print('Time not monotonic')
        time_diff_neg = np.nonzero(time_diff < 0)
        ind_diff_neg = time_diff_neg[0][0] + 1
        print('data file type, index of non-monotonic point = ', nn_mono, ind_diff_neg)
        err_file=open(errorfile,'a')    #Error log text file
        err_file.write(oformaterr1 % (datestring))
        err_file.write("%s" % " ")
        err_file.write(oformaterr2 % (nn_mono))
        err_file.write("%s" % " ")
        err_file.write(oformaterr3 % (ind_diff_neg))
        err_file.write("%s" % " ")
        err_file.write(oformaterr3 % (nvals_mono))
        err_file.write("%s" % " ")
        err_file.write(oformatmid1 % (time_mono[ind_diff_neg])) #Seconds
        err_file.write("%s" % "\n")
        err_file.close()

        #Put details into midnight file, to be read when processing other days
        #Only do this for co- and cross-polar files, not VAD
        out_range = nvals_mono - ind_diff_neg   #Number of points to write out
        midnight_file = midnight_file_root + datestring + '_' + str(nn_mono) + '.txt'
        #midnight_file = midnight_file_root + datestring + '_' + str(nn_mono) + '_v2.txt'
        mid_file = open(midnight_file,'w')
        for ka in range(out_range):
            mid_file.write(oformaterr1 % (datestring))
            mid_file.write("%s" % " ")
            mid_file.write(oformaterr2 % (nn_mono))
            mid_file.write("%s" % " ")
            mid_file.write(oformaterr3 % (ngates))
            mid_file.write("%s" % " ")
            mid_file.write(oformatmid1 % (time_mono[ind_diff_neg+ka]))
            mid_file.write("%s" % " ")
            mid_file.write(oformatmid1 % (elev_mono[ind_diff_neg+ka]))
            mid_file.write("%s" % "\n")
            for kb in range(ngates):    #Or is it ngates+1?
                mid_file.write(oformatmid1 % (vals_mono[(ind_diff_neg+ka),kb,0]))
                mid_file.write("%s" % " ")
                mid_file.write(oformatmid1 % (vals_mono[(ind_diff_neg+ka),kb,1]))
                mid_file.write("%s" % " ")
                mid_file.write(oformatmid1 % (vals_mono[(ind_diff_neg+ka),kb,2]))
                #If there are cleaned files, shouldn't be necessary now
                #if infile_type == 2:
                    #mid_file.write("%s" % " ")
                    #mid_file.write(oformatmid1 % (vals_cleaned[(ind_diff_neg+ka),kb,0]))
                    #mid_file.write("%s" % " ")
                    #mid_file.write(oformatmid1 % (vals_cleaned[(ind_diff_neg+ka),kb,1]))
                    #mid_file.write("%s" % " ")
                    #mid_file.write(oformatmid1 % (vals_cleaned[(ind_diff_neg+ka),kb,2]))
                mid_file.write("%s" % "\n")
        mid_file.close()

    else:       #Data all monotonic
        ind_diff_neg = -1

    return ind_diff_neg






# ------------------------------------------------------------------------
# Define the main netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_cl51(nday,cal_factor):
#Produces a netCDF file which meets AMF-CF requirements but can't produce a .png file using chilncplot from it 


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]
    #print(datevals)

    dat              = datetime.datetime(year,month,day,0,0,0)
    test_mon         = calendar.monthrange(year,month)
    #Test if it's the last day of the month - test_mon[1] has number of days in that month
    #if day < test_mon[1]:
        #next_dat     = datetime.datetime(year,month,(day+1),0,0,0)
    #else:
        #next_dat     = datetime.datetime(year,(month+1),1,0,0,0)
    start_of_day     = date2num(dat)
    daystr           = dat.strftime("%Y%m%d_")
    dirstr           = dat.strftime("%Y") + '/' + dat.strftime("%Y%m") + '/' + dat.strftime("%Y%m%d") + '/'
    #datestring       = dat.strftime("%Y%m%d")
    start_of_day_str = dat.strftime("%Y-%m-%d %H:%M")
    dat_str          = dat.strftime("%Y-%m-%dT%H:%M:%S") 
    #next_dat_str     = next_dat.strftime("%Y-%m-%dT%H:%M:%S")
    inpath           ="/data/vaisala_ceilometer/mirror_cl51sky_local_ctview/"
    infile_base      = "C" + datestring[3:]
    path_out = "/data/amof-netCDF/ncas-ceilometer-5/"
    outfile=path_out+datestring+"details.txt"
    errorfile = path_out + "cl51_errors.txt"
    oformaterr1="%s"
    oformaterr2="%7.3f %2i"

    # -------------------------
    # Find files for that day 
    # -------------------------
    files = os.listdir(inpath)
    expr    = re.compile(infile_base+'[0-2][0-9]\.DAT')
    infiles = [elem for elem in files if expr.search(elem)]
    print(infiles)
    infiles.sort()
    nfiles = len(infiles)
    if nfiles == 0:
        print('No files to process')
        #sys.exit()

    # ----------------------
    # Initialize data arrays
    # ----------------------
    n_index    = 5	#1-3 cloud bases, measured vertical vis (4), highest measured signal (4)
    n_sky_cond = 5	#5 cloud layers reported from sky condition algorithm
    fixed_az   = 280.0
    h_geoid    = 131.0	#Change - this is above msl
    timesecs   = np.zeros((20000))
    elev_angle = np.zeros((20000))
    az_angle   = np.zeros((20000))+fixed_az
    #inst_alarm = np.zeros((20000))
    braw       = np.zeros((20000,2000))
    nbase      = np.zeros((20000,n_index))		#As well as no. bases, is a diagnostic for vert. vis. details
    cbase      = np.zeros((20000,n_index))	#3 possible cloud base heights, plus vert. vis. if cloud status = 4. 2D to save as qc_flag
    cbase_write= np.zeros((20000,n_index))	#As cbase, but h_geoid added as it's read in, to avoid adding h_geoid to missing values
    cloudcov   = np.zeros((20000,n_sky_cond))	#Coverage of up to 5 cloud layes (oktas)
    cbcov      = np.zeros((20000,n_sky_cond))	#5 possible cloud base heights from sky condition algorithm
    qc_flag_cov= np.zeros((20000,n_sky_cond))	#Quality flag for cloud coverage
    inst       = np.zeros((20000,15))	#Instrument parameters
    #tempbgvar  = np.zeros((12))

    # -------------------------
    # Define various parameters
    # -------------------------
    hbin        = 10.0                  # Bin height (m), is in header information too
    snr_threshold = 3.0
    param = np.zeros((11))
    param[0] = 51.1445	#Latitude
    param[1] = 358.563	#Longitude
    param[2] = 84.0	#Height of site above geoid
    param[3] = 10000.0	#Valid max. altitude
    param[4] = 910.0	#Laser wavelength (nm)
    param[5] = 0.000003 #Pulse energy (J)
    param[6] = 6500     # PRF 
    param[7] = 0.148     # Lens diameter (m)
    param[8] = 0.014  # Transmit beam divergence (deg)
    param[9] = 1.1E-07  # Pulse length (s)
    param[10] = 3.0E+06  # Bandwidth of measurements at 3dB level (s-1)
    lyr_index = np.arange(n_index)	#Do we need this, or is there a way we can avoid including it as a variable?
    cov_index = np.arange(n_sky_cond)	#Do we need this, or is there a way we can avoid including it as a variable?
    assoc_file_string = "n/a" 
    missing_value = -1.0E+20
    missing_cloudcov_value = -127	#Not -999 as it's saved to netCDF file as a byte value. Avoiding 255 for now as it comes out as -1b 

    cbase[:,:] = missing_value	#Set default as missing as number of cloud base heights reported varies (0 - 3)
    cbase_write[:,:] = missing_value	#Value with h_geoid added for writing out, not processing
    cloudcov[:,:] = missing_cloudcov_value
    cbcov[:,:] = missing_value


    m = 0	#Data lines in all files from 1 day
    #f = open(infile, 'r')
    for nf in range(nfiles):	#Indent from here
        f = open(inpath+infiles[nf], 'r')
        #f = open(inpath+infiles[0], 'r')

        #First 2 lines appear 
        for z in range(2):
            line = f.readline()

        while True:
            err_sk  = 0                     #Flag for error in sky condition algorithm
            line = f.readline()	#Contains date and time
            #if not line: break
            if len(line) == 21:	#Date/time line, not "file closed" line at end 
                tempdata = line.split()	#Splits at space into date and time strings
                fdata=tempdata[1]
                timesecs[m] = 3600.0*float(fdata[0:2]) + 60.0*float(fdata[3:5]) + float(fdata[6:8])
            #else: break
            else:	#May be file closed/file created line
                for z in range(2):
                    line = f.readline()
                line = f.readline()	#
                if len(line) == 21:	#Date/time line, not "file closed" line at end 
               	    tempdata = line.split()	#Splits at space into date and time strings
               	    fdata=tempdata[1]
               	    timesecs[m] = 3600.0*float(fdata[0:2]) + 60.0*float(fdata[3:5]) + float(fdata[6:8])

            line = f.readline()	#No useful information
            if len(line) == 0: break
            line = f.readline()	#Alarms, warnings, options selected
            char1 = line[0:1]
            if char1.find('/') < 0:	#Should we look for /A or A as 1st character too?	
                nb = int(line[0])
            else:
                nb = 6	#Data quality suspect
                #Do I want to read cloud base heights if error? I think not, set as missing
                err_type = 1
                cbase[m,:] = missing_value 
                cbase_write[m,:] = missing_value
                print('Alarm condition found') 
                err_file=open(errorfile,'a')	#Error log text file
                err_file.write(oformaterr1 % (datestring))
                err_file.write("%s" % " ")
                err_file.write(oformaterr2 % (timesecs[m]/3600.0,err_type))
                err_file.write("%s" % "\n")
                err_file.close()

            tempdata = line.split()	#Haven't recorded alarms from this line. Choice of m/ft is implicit in "range" units.
            #print(line)
            #print(tempdata)
            if nb > 0 and nb < 4:
                for k in range(nb):
                    cbase[m,k]=int(tempdata[k+1])
                    cbase_write[m,k] = cbase[m,k] + h_geoid
            if nb == 4:
                cbase[m,3]=int(tempdata[1])	#If Cloud diag. = 4, vert. vis is given in place of 1st cloud base height
                cbase[m,4]=int(tempdata[2])	#If Cloud diag. = 4, highest detected signal is given in place of 2nd cloud base height
                cbase_write[m,3] = cbase[m,3] + h_geoid
                cbase_write[m,4] = cbase[m,4] + h_geoid
            nbase[m,:]=nb

            alarm_str=tempdata[4]
            inst[m,14] = int(alarm_str,16)
            #It will be hard to decode alarm string to get every alarm out. Can we just report it as now with ct75k? As a quality flag?
            #Check height units are metres
            if int(alarm_str[10:11]) != 8:
                print("Height units set to feet!!")
                #Should this cause a break function?
    
            line = f.readline()	#Cloud coverage
            tempdata = line.split()	#Leading spaces ignored, all strings at this stage
            if int(tempdata[0]) >= 0 and int(tempdata[0]) < 9:	#Valid values for no. of oktas coverage
                for ck in range(5):
                    if int(tempdata[(2*ck)]) > 0:	#& cbase[m,ck] > 0	#Possibly only write coverage if have a cbase value?
                        cloudcov[m,ck]=int(tempdata[(2*ck)])	#Only writes cloud coverage data to variable if it's > 0. Correct?
                        cbcov[m,ck]=10*int(tempdata[(2*ck+1)]) + h_geoid	#cbcov doesn't get used in data processing
                qc_flag_cov[m,:] = 1
            else:   #i.e. tempdata[0] = 99 or -1
                if int(tempdata[0]) == 99:
                    qc_flag_cov[m,:] = 2
                    print(m, 'Sky condition algorithm has insufficient data')
                    err_type = 2
                    err_file=open(errorfile,'a')	#Error log text file
                    err_file.write(oformaterr1 % (datestring))
                    err_file.write("%s" % " ")
                    err_file.write(oformaterr2 % (timesecs[m]/3600.0,err_type))
                    err_file.write("%s" % "\n")
               	    err_file.close()
                if int(tempdata[0]) == -1:
                    qc_flag_cov[m,:] = 3
                    print(m, 'Sky condition algorithm not operating')
                    err_type = 3
                    err_file=open(errorfile,'a')	#Error log text file
                    err_file.write(oformaterr1 % (datestring))
                    err_file.write("%s" % " ")
                    err_file.write(oformaterr2 % (timesecs[m]/3600.0,err_type))
                    err_file.write("%s" % "\n")
                    err_file.close()
		     
                err_sk = 1
                #nbase[m,:] = 6        #Error in measurement
 
            line = f.readline()	#Instrument data
            tempdata = line.split()
            ngates=int(tempdata[2])
            for ni in range(8):
                inst[m,ni] = float(tempdata[ni])	#Scaling, res., samples, pulse energy, temp, transmission, tilt, bg
            instr = tempdata[8]
            if instr[0] == 'L':
                inst[m,8] = 1
            else:
                inst[m,8] = 2
            inst[m,9] = 1024.0*float(instr[1:5])
            if instr[5] == 'H':
                inst[m,10] = 1
            else:
                inst[m,10] = 2	#Was zero, changed to match AMF quality flag requirements
            if instr[6] == 'N':
                inst[m,11] = 1
            else:
                inst[m,11] = 2
            inst[m,12] = float(instr[7:])
            inst[m,13] = float(tempdata[9])

            line = f.readline()	#Backscattering profile
            #if err_sk == 0: #No error in sky condition algorithm
            for n in range(ngates):
                tempval=int(line[(5*n):(5*n)+5],16)
                # Invert data where values are ~2^20, or 1048576
                if tempval >= 1.0E+06:
                    tempval = tempval - 1048576
                braw[m,n]=float(tempval)
            if np.amax(braw[m,:]) <= 1.0:	#Still in raw units here, this is a small backscattering coefficient, for if no profile values are recorded
                braw[m,:] = missing_value
                #nbase[m,:] = 6
                print('Profile missing')
                err_type = 3
                err_file=open(errorfile,'a')	#Error log text file
                err_file.write(oformaterr1 % (datestring))
                err_file.write("%s" % " ")
                err_file.write(oformaterr2 % (timesecs[m]/3600.0,err_type))
                err_file.write("%s" % "\n")
               	err_file.close()


            line = f.readline()	#No useful information
            line = f.readline()	#Empty line
            m += 1	

        f.close()
        #Indent to here, for reading each file from day

    if m > 0:	#1 or more data lines read

        end_m = m
        #Look for a last time stamp that is lower than the previous (occasionally last time is 00:00:00 from next day
        if timesecs[m-2] > timesecs[m-1]:
            end_m = m-1
        #Reduce array sizes
        timesecs = timesecs[0:end_m]
        braw = braw[0:end_m,0:ngates]
        cbase = cbase[0:end_m,:]
        cbase_write = cbase_write[0:end_m,:]
        cbase_positive = cbase_write[np.where(cbase_write >= 0.0)]	#Used for finding minimum cbase value in data
        if len(cbase_positive) > 0:
            min_cbase_positive = np.amin(cbase_positive)
            max_cbase_positive = np.amax(cbase_positive)
        else:
            min_cbase_positive = missing_value
            max_cbase_positive = missing_value
        nbase = nbase[0:end_m,:]
        nbase_write = nbase + 1	#0 isn't used in quality flag, so no cloud bases = 1, 1 cloud base = 2, etc. Ugly!
        cloudcov = cloudcov[0:end_m,:]
        cloudcov_positive = cloudcov[np.where(cloudcov >= 0.0)]	#Used for finding minimum cloud coverage value in data
        if len(cloudcov_positive) > 0:
            min_cloudcov_positive = np.amin(cloudcov_positive)
            max_cloudcov_positive = np.amax(cloudcov_positive)
        else:
            min_cloudcov_positive = missing_cloudcov_value
            max_cloudcov_positive = missing_cloudcov_value
        cbcov = cbcov[0:end_m,:]
        cbcov_positive = cbcov[np.where(cbcov >= 0.0)]	#Used for finding minimum cb value in data
        if len(cbcov_positive) > 0:
            min_cbcov_positive = np.amin(cbcov_positive)
            max_cbcov_positive = np.amax(cbcov_positive)
        else:
            min_cbcov_positive = missing_value
            max_cbcov_positive = missing_value
        qc_flag_cov = qc_flag_cov[0:end_m,:]
        inst = inst[0:end_m,:]
        az_angle = az_angle[0:end_m]
        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        elev_angle=inst[:,6]	#from zenith
        gain_flag = np.broadcast_to(inst[:,10],(ngates,end_m))
        gain_flag = np.transpose(gain_flag)


        max_braw=np.amax(braw)
        print('Maximum braw = ',max_braw)
        braw = 1.0E-08*braw	#Scaling to m-1 sr-1 from units in raw file.Doesn't fit nc variable if transpose.
        #braw_real = braw[np.where(braw > missing_value)]	#Used for finding minimum braw value in data
        braw_real = braw[np.where(braw >= 0.0)]	#Used for finding minimum braw value in data

        # ---------------------
        # Setup height variable
        # ---------------------
        h    = np.zeros((ngates))
        h[:] = range(ngates)
        h    = h + 0.5	#To refer height to centre of bin, not bottom
        h *= hbin * 0.001 

        #--------------------------------------------
        # Calculating background and removing speckle
        #--------------------------------------------


        #Remove range^2 correction
        braw_nor2=braw/(h*h)
        #Possibly need to do b/g subtraction again here?

        #Calculate background for each profile from average of last 100 bins, approx. 1 km
        #bgvar=np.var(braw_nor2[:,401:500],axis=1)	#Try working out variance at lower altitude, looks noisier when cloud is low   
        bgvar=np.var(braw_nor2[:,(ngates-99):ngates],axis=1)	#Simplest way to work out bgvar, but not great for low cloud

        aa_arr=np.arange(15)
        keep_sizes=np.zeros(end_m)
        aa_min_arr = ngates - (((aa_arr+1)*100)-1)
        for jj in range(end_m):
            #if (cbase[jj,0] > 0 and cbase[jj,0] <= 1000 and nbase[jj,0] < 4 and np.amax(cbase[jj,:]) <= 3000) or nbase[jj,0] == 4 or inst[jj,13] > 300:	#Was cbase[jj,0] < 500	#v1
            if (cbase[jj,0] > 0 and nbase[jj,0] < 4 and np.amax(cbase[jj,:]) <= 3000) or nbase[jj,0] == 4 or inst[jj,13] > 300:	#v2
                #Could change to (nbase in range 1-4 and max(cbase) < 3000) or sumbackscatter >300 
                max_cloud = np.amax(cbase[jj,:])+100	#Add 100m to be sure we're above cloud - started as 500m 
                max_cloud_index=np.argmin(np.absolute((1000*h)-max_cloud))
                pos_aa_min_arr = aa_min_arr[np.where((aa_min_arr-max_cloud_index) > 0)]
                size_arr=np.size(pos_aa_min_arr)
                keep_sizes[jj] = size_arr 
                tempbgvar=np.zeros(size_arr)	#Define for each profile as size can change 
                for aa in range(size_arr):	#Was 12, before starting to work out maximum index
                    low_ind = pos_aa_min_arr[aa]
                    high_ind = low_ind+99
                    tempbgvar[aa]=np.var(braw_nor2[jj,low_ind:high_ind])
                bgvar[jj] = np.amax(tempbgvar)
   
        snr = np.divide(np.transpose(braw_nor2),np.sqrt(bgvar))	#Transposed as otherwise division is rejected due to relative array shapes
        snr = np.transpose(snr)

        bs_flag = np.zeros((end_m,ngates))	#Was zeros, now 1 if poor, 0 if good
        n_masked = 0
        for kk in range(end_m):
            for ll in range(ngates):
                #if snr[kk,ll] >= snr_threshold and inst[kk,10] == 1 and np.sqrt(bgvar[kk]) > 3.0E-08:	#Data not included if low gain
                #if snr[kk,ll] >= snr_threshold and inst[kk,10] == 1:	#Data not included if low gain
                if snr[kk,ll] >= snr_threshold:	#Process for all data, high or low gain
                    if inst[kk,10] == 1: 
                        bs_flag[kk,ll] = 1	#Good data high gain high snr
                    else:
                        bs_flag[kk,ll] = 2	#Low gain high snr
                    #Look at surrounding points to check if this point is speckle
                    if kk > 0 and kk < (end_m-1) and ll > 0 and ll < (ngates-1):
                        array_to_check = np.asarray([snr[kk-1,ll], snr[kk+1,ll], snr[kk,ll-1], snr[kk,ll+1]])
                        #if np.amax(array_to_check) < snr_threshold:
                            #bs_flag[kk,ll] = 6
                            #n_masked = n_masked + 1
                        yes_no_check =np.where(array_to_check < snr_threshold, 1, 0)	#Check whether each of 4 adjacent points is below SNR threshold
                        if round(np.sum(yes_no_check)) >= 3:	#If at least 3 points are below SNR threshold
                            bs_flag[kk,ll] = 6
                            n_masked = n_masked + 1
                             
                else:				#low snr
                    if inst[kk,10] == 1: 
                        bs_flag[kk,ll] = 3	#High gain low snr
                    else:
                        bs_flag[kk,ll] = 4	#Low gain low snr
		
                    #bgvar filtering needs to be in different place - takes out whole profile for that time here
                    #The low SD would mean a higher value of snr
                if braw[kk,ll] < 0.0:
                    bs_flag[kk,ll] = 5	#Signal less than 0 - usually instrument noise
            if np.amax(braw[kk,:]) == missing_value:
                bs_flag[kk,:] = 7	#Profile missing
        print('Number of speckle points qc-flagged = ', n_masked)


        #--------------------------------------------
        # Adding height of geoid to cloud coverage 
        #--------------------------------------------

        h *= 1000*cos(0.01745*inst[0,6])	#Could make it mean of inst[:,6] through day
        h = h + h_geoid	#To be above WGS84 geoid
        #cbase = h_geoid + cbase
        #cbcov = h_geoid + cbcov

        #--------------------------------------------
        # Combining gain and backscattering flags 
        #--------------------------------------------


        #For backscatter file, gain and backscattering flags can be added
        #0 = low gain and low SNR, 1 = high gain, low SNR, 2 = high gain, high SNR.
        #If gain is low, bs_flag is automatically set to zero earlier
        #combined_flag = gain_flag + bs_flag
	

        #--------------------------------------------
        # Sorting out day/time information 
        #--------------------------------------------

        tt = (year,month,day,0,0,0,0,0,0)
        tt_upper = (year,month,(day+1),0,0,0,0,0,0)
        day_start_epoch = time.mktime(tt)
        #t_valid_min = day_start_epoch
        #t_valid_max = time.mktime(tt_upper) 
        yr_arr=np.zeros(end_m)+year
        mn_arr=np.zeros(end_m)+month
        dy_arr=np.zeros(end_m)+day
        dyyr_arr=np.zeros(end_m)
        epoch_timesecs=np.zeros(end_m)
        hr_arr=np.zeros(end_m,dtype=np.int32)
        mi_arr=np.zeros(end_m,dtype=np.int32)
        sc_arr=np.zeros(end_m,dtype=np.float32)
        epoch_timesecs = day_start_epoch + timesecs
        #In cases where parts of previous or next day are found in file, it would be best to derive the above from the data timestamp
        for mm in range(end_m):
            #epoch_timesecs[mm] = day_start_epoch + timesecs[mm]
            tt = time.gmtime(epoch_timesecs[mm])
            hr_arr[mm] = np.int32(tt[3])
            mi_arr[mm] = np.int32(tt[4])
            #sc_arr[mm] = np.float32(tt[5])
            sc_arr[mm] = timesecs[mm]-3600*hr_arr[mm]-60*mi_arr[mm]	#Keeps digits after decimal place, unlike tt[5] 
            dyyr_arr[mm] = np.float32(tt[7]+timesecs[mm]/86400.0)
            #Also deal with profiles with instrument alarms in this loop
            if nbase[mm,0] == 6:	#Data quality suspect
                bs_flag[mm,:] = 5
                nbase_write[mm,:] = 6

        #------------Setting up YAML data object----------------------
        template_file_path1 = "/home/chilbolton_software/python/ncas_python/halo/ncas-ceilometer-5-template_cbh.yaml"	#YAML template
        template_file_path2 = "/home/chilbolton_software/python/ncas_python/halo/ncas-ceilometer-5-template_cov.yaml"	#YAML template
        template_file_path3 = "/home/chilbolton_software/python/ncas_python/halo/ncas-ceilometer-5-template_rawbs.yaml"	#YAML template
        template_file_path4 = "/home/chilbolton_software/python/ncas_python/halo/ncas-ceilometer-5-template_calbs.yaml"	#YAML template
        tfp = [template_file_path1, template_file_path2, template_file_path3, template_file_path4]
        #tfp = [template_file_path2, template_file_path1]
        data_name = ["_cloud-base_v1.0", "_cloud-coverage_v1.0", "_aerosol-backscatter_standard_v1.0", "_aerosol-backscatter_advanced_v1.0"]
        cfarr_head = 'ncas-ceilometer-5_cao_'

        #--------------------------------------------
        # History information 
        #--------------------------------------------

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        conversion_time = datetime.datetime.utcnow()
        history_string = dat_str + " - data collection commenced using Chilbolton CL51 ceilometer s/n K0150001.\n"
        history_string = history_string + conversion_time.strftime("%Y-%m-%dT%H:%M:%S") + " - converted to netCDF from raw data by " + user_id +" on " + computer_id + "."


        #Only write out calibrated files if cal_factor > 0. It's set to zero if you don't want to write calibrated files
        if cal_factor > 0.0:
            total_files = 4
        else:
            total_files = 3

        for n_outfile in range(total_files):

            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[n_outfile])
            #data_object_id = creator.load_a_template(template_file_path)
            print('Data object ID = ', data_object_id)
            #handler.show_requirements_for_template(data_object_id)

            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------

            out_file   = os.path.join(path_out,cfarr_head+datestring+data_name[n_outfile]+'.nc')
            print('Output NetCDF file ' + out_file)

            #lengths_of_dimensions depends on which file you're writing, as layer_index is used differently in cbh and coverage

            if n_outfile == 0:	#CBH
                lengths_of_dimensions = {"time": end_m, "layer_index": n_index}

            if n_outfile == 1:	#Coverage
                lengths_of_dimensions = {"time": end_m, "layer_index": n_sky_cond}

            if n_outfile == 2:	#Uncalibrated beta values
                lengths_of_dimensions = {"time": end_m, "altitude": ngates}
                assoc_file_string = "Instrument is at ground level.\nThe time variable gives the time at the end of the acquisition interval."

            if n_outfile == 3:	#Calibrated data
                lengths_of_dimensions = {"time": end_m, "altitude": ngates}
                braw_real = cal_factor * braw_real
                braw = cal_factor * braw
                braw = np.where(braw > cal_factor * missing_value, braw, missing_value)		    
                assoc_file_string = "Instrument is at ground level.\nThe time variable gives the time at the end of the acquisition interval.\n" + "Associated file: " + cfarr_head + datestring + data_name[2] + ".nc \nwhich contains attenuated_aerosol_backscatter_coefficient as recorded by the ceilometer."
                history_string = history_string + "\nand attenuated_aerosol_backscatter_coefficient calibrated by multiplying by a factor of " + str(cal_factor) + "."

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],int(sc_arr[0])), "end_time": datetime.datetime(year,month,day,hr_arr[end_m-1],mi_arr[end_m-1],int(sc_arr[end_m-1])), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec": np.amax(sc_arr), "h_min": np.amin(h), "h_max": np.amax(h), "zen_min": np.amin(inst[:,6]), "zen_max": np.amax(inst[:,6]), "az_min": np.amin(az_angle), "az_max": np.amax(az_angle), "min_pulses": np.amin(inst[:,9]), "max_pulses": np.amax(inst[:,9]), "min_p_scaling": np.amin(inst[:,0]), "max_p_scaling": np.amax(inst[:,0]), "min_en_scaling": np.amin(inst[:,3]), "max_en_scaling": np.amax(inst[:,3]), "min_laser_T": (np.amin(inst[:,4])+273.0), "max_laser_T": (np.amax(inst[:,4])+273.0), "min_trans": np.amin(inst[:,5]), "max_trans": np.amax(inst[:,5]), "min_bg": np.amin(inst[:,7]), "max_bg": np.amax(inst[:,7]), "min_samp": (np.amin(inst[:,12])*1000000), "max_samp": (np.amax(inst[:,12])*1000000), "min_alarm": np.amin(inst[:,14]), "max_alarm": np.amax(inst[:,14]), "min_bs_sum": (np.amin(inst[:,13])*0.0001), "max_bs_sum": (np.amax(inst[:,13])*0.0001), "min_bs": np.amin(braw_real), "max_bs": np.amax(braw_real), "min_cb": min_cbase_positive, "max_cb": max_cbase_positive, "min_cbcov": min_cbcov_positive, "max_cbcov": max_cbcov_positive, "min_cloudcov": min_cloudcov_positive, "max_cloudcov": max_cloudcov_positive, "comment_str": assoc_file_string, "history_str": history_string, "uname": user_id, "codename": __file__, "machinename": computer_id}


            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]
            data_object["variables"]["sensor_zenith_angle"]["values"][:] = inst[:,6]
            data_object["variables"]["sensor_azimuth_angle"]["values"][:] = az_angle[:] 
            data_object["variables"]["profile_pulses"]["values"][:] = inst[:,9]
            data_object["variables"]["profile_scaling"]["values"][:] = inst[:,0]
            data_object["variables"]["laser_pulse_energy"]["values"][:] = inst[:,3]
            data_object["variables"]["laser_temperature"]["values"][:] = inst[:,4]+273.0
            data_object["variables"]["window_transmittance"]["values"][:] = inst[:,5]
            data_object["variables"]["background_light"]["values"][:] = inst[:,7]
            #data_object["variables"]["sampling_frequency"]["values"][:] = inst[:,12]*1000000 
            #data_object["variables"]["instrument_alarm"]["values"][:] = inst[:,14]
            data_object["variables"]["backscatter_sum"]["values"][:] = inst[:,13]*0.0001

            if n_outfile == 0:
                #Cloud base height
               	#data_object["variables"]["layer_index"]["values"][:] = lyr_index[:]
               	data_object["variables"]["cloud_base_altitude"]["values"][:] = cbase_write[:]	
                data_object["variables"]["qc_flag"]["values"][:] = nbase_write[:]	
		
            if n_outfile == 2 or n_outfile == 3:
                #1=Raw aerosol backscatter, 2=Calibrated aerosol backscatter, may need to be separate eventually to use different yaml files
               	data_object["variables"]["altitude"]["values"][:] = h[:]
               	data_object["variables"]["attenuated_aerosol_backscatter_coefficient"]["values"][:] = braw[:]	#Why not [:,:]?	
               	data_object["variables"]["qc_flag"]["values"][:] = bs_flag[:]	
                #data_object["variables"]["qc_flag_gain"]["values"][:] = inst[:,10]	

            if n_outfile == 1:
                #Cloud coverage from sky condition algorithm 
                data_object["variables"]["cloud_base_altitude"]["values"][:] = cbcov[:]	
                data_object["variables"]["cloud_coverage"]["values"][:] = cloudcov[:]	
               	data_object["variables"]["qc_flag"]["values"][:] = qc_flag_cov[:]	
                #data_object["variables"]["sky_condition_index"]["values"][:] = cov_index[:]
                #data_object["variables"]["layer_index"]["values"][:] = cov_index[:]	#How NCAS want it, but a pain having dimension same name in 2 files
        
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print('Exit code = ', exit_code)

            print('Changing permissions of file ',out_file)
            print(' ')
            #oscommand = "chmod g+w " + out_file
            oscommand = "chmod 775 " + out_file
            os.system(oscommand)

            # Plotting data
            # Can't be done until file is written
            if n_outfile == 2:
                lidar_plot.lidar_plotting(nday, 1, 1)   #Plot raw data
            if n_outfile == 3:
                lidar_plot.lidar_plotting(nday, 1, 2)   #Plot calibrated data
            if n_outfile == 0:
                lidar_plot.lidar_plotting_time_series(nday, 1, 3)
            if n_outfile == 1:
                lidar_plot.lidar_plotting_time_series(nday, 1, 4)

        #-------------Text file for writing time series variables to separate file-----------------------------
        ofile=open(outfile,'w')	#Details text file
        oformat1="%7.3f %10.3e %2i %4i %2i"
        oformat2="%7.3f %10.3e %2i %10.3e %4i %4i %4i %4i %4i %2i %2i"
        for k in range(0,end_m):
            #ofile.write(oformat1 % (timesecs[k]/3600.0,np.sqrt(bgvar[k]),inst[k,10],cbase[k,0],nbase[k]))
            ofile.write(oformat2 % (timesecs[k]/3600.0,np.sqrt(bgvar[k]),inst[k,10],inst[k,13],cbase[k,0],cbase[k,1],cbase[k,2],cbase[k,3],cbase[k,4],nbase[k,0],keep_sizes[k]))
            ofile.write("%s" % "\n")
        ofile.close()


# ------------------------------------------------------------------------
# Define the general netcdf generation function, common to all sensors
# ------------------------------------------------------------------------
#def generate_netcdf_common(nday,year,month,day,datestring):
def generate_netcdf_common(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------
    ndat             = num2date(nday)
    year             = ndat.year
    month            = ndat.month
    day              = ndat.day
    dat              = datetime.datetime(year,month,day,0,0,0)
    nowdat           = datetime.datetime.now()
    nowstring        = nowdat.strftime("%Y-%m-%d %H:%M:%S")
    datestring       = dat.strftime("%Y%m%d")
    start_of_day_str = dat.strftime("%Y-%m-%d %H:%M")

    return year,month,day,datestring,start_of_day_str,nowstring


# ------------------------------------------------------------------------
# Define a function to broadcast a 2D array to 3D
# ------------------------------------------------------------------------
def my_broadcast(data_arr,no_t_av,ngates,nangles):

    d_arr_int1=np.broadcast_to(data_arr,(nangles, no_t_av, ngates))
    d_arr_int2=np.swapaxes(d_arr_int1,0,2)
    d_arr_out=np.swapaxes(d_arr_int2,0,1)
    return d_arr_out
