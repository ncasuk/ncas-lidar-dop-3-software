function ret = ceil_halo_calibrate(yyyymmdd_start,yyyymmdd_end,thres)
%Finds suitable data for calibrating the ceilometer and Halo lidar and calculates a calibration factor 
%Judith Agnew 20130404 
%Modified to put more of code in a function and add new Halo lidar
%Judith Jeffery 20180110
%Modified to read Halo2 files from Ewan's script, generated automatically

if (nargin < 3)
        disp('usage: ceil_halo_calibrate(yyyymmdd at start, yyyymmdd at finish, cloud threshold (1 or 2 for 1e-04 or 2e-04)');
        return;
end

if not(ischar(yyyymmdd_start))
        yyyymmdd_start = num2str(yyyymmdd_start);
end
if not(ischar(yyyymmdd_end))
        yyyymmdd_end = num2str(yyyymmdd_end);
end
if not(isnumeric(thres))

        thres = str2num(thres);
end


load('/home/jla/matlab/hoganjet.dat');
hoganjet=hoganjet'/255.0;
colormap(hoganjet);


year1=str2num(yyyymmdd_start(1:4));
month1=str2num(yyyymmdd_start(5:6));
day1=str2num(yyyymmdd_start(7:8));
year2=str2num(yyyymmdd_end(1:4));
month2=str2num(yyyymmdd_end(5:6));
day2=str2num(yyyymmdd_end(7:8));

%Default values - needed in main code as used in filenames
%cloud_thres=2e-04;	%threshold for a liquid cloud, default 2e-04, started at 8e-05
%bl_thres=1.0e-05;	%Drizzle/aerosol threshold for accepting point, start at 1.0e-05
%Current values
cloud_thres=zeros(1,4)
bl_thres=zeros(1,4)
cloud_thres(1:4) = thres * 1.0e-04;
bl_thres(1:4)=1.0e-05;
peakstr=num2str(cloud_thres(1),'%1.0e');
blstr=num2str(bl_thres(1),'%1.0e');
%bin_step=0.1;
%bin=0.05:bin_step:10;	%bins for calibration factor values
bin_step=0.05;
bin=0.025:bin_step:5;	%bins for calibration factor values
pdfcalc=zeros(2,length(bin),13);	%pdf array 2 x instruments x bins x months

%File for statistics
%Shows number of profiles in day, number with biggest signal in profile 0.5-2.0km
%sum of number of profiles in an hour with>=95% of peaks above threshold, number with background below cloud less than threshold
out_stats_file=['/data/range/ceil_halo_cal/ceil_halo_stats_',yyyymmdd_start,'_',yyyymmdd_end,'_',peakstr,'_',blstr,'.txt'];
fid_stats=fopen(out_stats_file,'w');

start_ts = datenum( year1, month1, day1 , 0, 0, 0);
finish_ts = datenum( year2, month2, day2, 0, 0, 0);

current_ts = start_ts;
for current_ts = start_ts:finish_ts
        a = datestr(current_ts, 30);
        now_date = a(1:8)
	num_month=str2num(a(5:6))
	m=1;	%Counter for integrated peaks
	fprintf(fid_stats,'%s\t',now_date);

%-----CT75k ceilometer file details, now defunct--------------------------------------------------
	%if current_ts >= datenum(2008, 9, 1, 0, 0, 0)
	if current_ts >= datenum(2006, 1, 1, 0, 0, 0)
		%Always the case since files have been recalibrated
		inpath=['/home/eoconnor/calibrated/ct75k/',now_date(1:4),'/']; 
		fname1=[inpath, now_date, '_chilbolton-ct75k.nc'];
	elseif current_ts <= datenum(2004,2,29,0,0,0)
		fname1=['/home/jla/ct75k/lidar-ct75k_chilbolton_',now_date, '.nc'];
	else
		fname1=['/home/jla/ct75k/cfarr-lidar-ct75k_chilbolton_',now_date, '.nc'];
	end
	%out_ceil_file=['/home/jla/lidar_cal/ceil_cal_',yyyymmdd,'.txt']; 
%-----Original mk1 HALO file details, now defunct--------------------------------------------------
	if current_ts >= datenum(2010, 9, 1, 0, 0, 0)
		inpath=['/home/halo/calibrated/halo-doppler-lidar/',now_date(1:4),'/'];	%should I use tilt data to closer match ceilometer path? 
		fname2=[inpath, now_date, '_chilbolton_halo-doppler-lidar.nc'];
	else
		inpath=['/home/jla/halo/']
		fname2=[inpath, 'cfarr-lidar-doppler-halo_chilbolton_',now_date,'.nc'];
	end
	if current_ts >= datenum(2014,1,16,0,0,0) & current_ts <= datenum(2014,4,15,0,0,0)
		inpath='/home/jla/lidar_cal/loan_halo/';
		fname2=[inpath,now_date,'_chilbolton_halo-doppler-lidar.nc'];
	end	
	%out_halo_file=['/home/jla/lidar_cal/halo_cal_',yyyymmdd,'.txt']; 
%-----CL51 ceilometer file details--------------------------------------------------
	inpath=['/home/eoconnor/uncalibrated/cl51/',now_date(1:4),'/'];
	fname3=[inpath,now_date,'_chilbolton_cl51.nc']
%-----Halo2 file details--------------------------------------------------
        %Change which lines are commented out here to swap between using files containing as-recorded data (uncleaned)
	%or files processed using the FMI toolbox (cleaned)
	%Cleaned files give a more consistent calibration value
	%inpath=['/data/lidar-doppler-halo/netcdf/uncalibrated/halo-doppler-lidar-118-co/',num2str(year1),'/'];	%Uncleaned Halo file directory
	inpath=['/data/lidar-doppler-halo/netcdf/calibrated/halo-doppler-lidar-118-co/',num2str(year1),'/'];	%Cleaned Halo file directory
	%fname4=[inpath,now_date,'_chilbolton_halo-doppler-lidar-118-co.nc']	%Uncleaned Halo filename
	fname4=[inpath,now_date,'_chilbolton_halo-doppler-lidar-118-stare-co.nc']	%Uncleaned Halo filename
	%fname=char(fname1,fname2,fname3,fname4);	%CT75K, Old Halo, CL51, Halo2
	fname=char(fname3,fname4);	%CL51, Halo2

	%Changed for working with loan Halo lidar as this period starts and ends mid-month
	%ceil_plot_file=['/data/range/ceil_halo_cal/ceil_cal_newceilcal_',yyyymmdd_start,'_',yyyymmdd_end,'.jpg'];

	halo_plot_file=['/data/range/ceil_halo_cal/ceil_halo_cal_',yyyymmdd_start,'_',yyyymmdd_end,'_',peakstr,'_',blstr,'.png'];
	out_text_file=['/data/range/ceil_halo_cal/ceil_halo_pdf_',yyyymmdd_start,'_',yyyymmdd_end,'_',peakstr,'_',blstr,'.txt'];

%------Test that all files exist-----------------------------------------

	fid3=fopen(deblank(fname(1,:)),'r')
	fid4=fopen(deblank(fname(2,:)),'r')
	t(1)=fid3;
	t(2)=fid4;
	if fid3 > 0
		fclose(fid3);
	end
	if fid4 > 0
		fclose(fid4);
	end

	val = calculate_calibration(fname,cloud_thres,bl_thres,t,fid_stats,num_month,m,bin_step,bin);

	fprintf(fid_stats,'\n');	%Outside t1, t2, t3>0 loop so that it happens even if a file doesn't exist

        current_ts = current_ts + datenum (0, 0, 1, 0, 0, 0);

	pdftemp = val.pdfout;	%2D pdf array from function

	max(max(max(pdfcalc)))

	pdfcalc(:,:,num_month) = pdfcalc(:,:,num_month) + pdftemp;
	sum(pdfcalc(1,:,num_month))
	sum(pdfcalc(2,:,num_month))
	%sum(pdfcalc(3,:,num_month))

end	%days loop
	fclose(fid_stats);

	fid=fopen(out_text_file,'w');
	
	for n=1:12
		fprintf(fid,'%d\t',n);
		for k=1:2
			for m=1:length(bin)
				fprintf(fid,'%d\t',pdfcalc(k,m,n));
			end
		end
	fprintf(fid,'\n');
	end

	fclose(fid);



	subplot(2,1,1);
	pcolor(1:13,bin,squeeze(pdfcalc(1,:,:))); shading flat
	title(['CL51 calibration plot ',yyyymmdd_start(1:4)]); xlabel('Month'); ylabel('Calibration factor'); xlim([1 13]); ylim([0 4]);
	subplot(2,1,2);
	pcolor(1:13,bin,squeeze(pdfcalc(2,:,:))); shading flat
	title(['Halo 2 calibration plot ',yyyymmdd_start(1:4)]); xlabel('Month'); ylabel('Calibration factor'); xlim([1 13]); ylim([0 4]);
	%[img,cmap]=capture;	%Older version of Matlab
	%imwrite(img,cmap,halo_plot_file,'Quality',85);	%Older version of Matlab
	saveas(gcf,halo_plot_file)

%--------------------------------------------------------------------------------------------------------
function val = calculate_calibration(fname,cloud_thres,bl_thres,t,fid_stats,num_month,m,bin_step,bin,pdf)
%--------------------------------------------------------------------------------------------------------
% Takes backscattering data and calculates calibration factor

pdf=zeros(2,length(bin));	%pdf array instruments * bins
sum_cloud_ok=0;	%For stats of which clouds are included

%Coefficients for eta*S
A=0.02704;
B=-0.3469;
C=1.681;
D=-4.142;
E=18.57;
h_etaS=20.0;

height_scale(1)=0.99756;
height_scale(2)=0.001;
bl_low_ind(1)=2;
bl_low_ind(2)=3;
bl_high_ind(1)=10;
bl_high_ind(2)=2;
width(1)=9;
width(2)=3;

for nn=1:2

	if t(nn) > 0

		nc=netcdf(deblank(fname(nn,:)),'nowrite');
		len=length(nc);	%check there are data in the file
		if len > 0
			time_temp=nc{'time'};
			ceil_time=time_temp(:);
			disp(['Number of time values = ',num2str(length(ceil_time))]);
			height_temp=nc{'range'};
			if length(ceil_time) > 0 & length(height_temp(:)) > 0	%Are arrays complete in file
			ceil_height=height_scale(nn)*height_temp(:);
			bscatt_temp=nc{'beta'};
			if nn == 1 	%Ceilometer files use scale factor
				scale_fac=nc{'beta'}.scale_factor(:);
				ceil_bscatt=scale_fac*bscatt_temp(:);
			else
				ceil_bscatt=bscatt_temp(:);	%number is just a temporary calibration factor
			end
			ceil_bscatt(2,1:20)
			size(ceil_bscatt)
			ceil_h_step=1000.0*(ceil_height(2)-ceil_height(1));	%Convert to m
			[max_ceil,ceil_index]=max(ceil_bscatt,[],2);	%ceil_index=height indicies of max. bscatt values
			c_heights_1=ceil_height(ceil_index);	%heights of max. bscatt values, still for all times
			%c_heights_1(1:20)
			fprintf(fid_stats,'%d\t',length(ceil_time));
			ceil_inrange=find(c_heights_1 >= 0.5 & c_heights_1 <= 2.0);
			length(ceil_inrange)
			fprintf(fid_stats,'%d\t',length(ceil_inrange));

%-----------------------Check hour by hour for enough points above cloud threshold---------------------------
%			Work from ceil_inrange subsets
			c_time_1=ceil_time(ceil_inrange);
			c_heights_1=c_heights_1(ceil_inrange);
			c_maxbscatt_1=max_ceil(ceil_inrange);
			
			for tlow=0:23
				[junk,low_ind]=min(abs(c_time_1-tlow));
				[junk,high_ind]=min(abs(c_time_1-tlow-1));
				[junk,low_ind_all]=min(abs(ceil_time-tlow));
				[junk,high_ind_all]=min(abs(ceil_time-tlow-1));
				high_ind=high_ind-1;	%go to 1 lower index to avoid overlap
				high_ind_all=high_ind_all-1;	%go to 1 lower index to avoid overlap
				n_indicies=high_ind_all-low_ind_all+1;
				cloud_ok=find(c_maxbscatt_1(low_ind:high_ind) >= cloud_thres(nn));
				ratio=length(cloud_ok)/n_indicies;
				if ratio < 0.95	%not enough of the hour above the cloud threshold
					c_maxbscatt_1(low_ind:high_ind)=-0.0001;	%Will search for individual points below the cloud threshold later and catch these
				end
				sum_cloud_ok=sum_cloud_ok+length(cloud_ok);
			end	%tlow for loop
			fprintf(fid_stats,'%d\t',sum_cloud_ok);
			sum_cloud_ok=0;	%Will be ready for next use

%---------------------- Look for few remaining values below cloud threshold-------------------------------

			cloud_ok=find(c_maxbscatt_1 >= cloud_thres(nn));
			fprintf(fid_stats,'%d\t',length(cloud_ok));
			c_time_1=c_time_1(cloud_ok);
			c_heights_1=c_heights_1(cloud_ok);
			c_maxbscatt_1=c_maxbscatt_1(cloud_ok);
			n_c_points=length(c_time_1);

			%New arrays for acceptable peaks for integration, max. size as previous arrays
			c_time_2=zeros(1,n_c_points);
			c_calfactor=zeros(1,n_c_points);		

%----------------------Loop through remaining points, checking for drizzle/high aerosol and integrating peaks----------

			for n=1:n_c_points
				%Find profile corresponding to each acceptable point, up to 4 points below the peak
				[junk,prof_num]=min(abs(ceil_time-c_time_1(n)));
				[junk,c_h_ind]=min(abs(ceil_height-c_heights_1(n)));
				ext_peak=ceil_bscatt(prof_num,(c_h_ind-width(nn)):(c_h_ind+width(nn)));
				peak_points=find(ext_peak >= 0.05*c_maxbscatt_1(n));
				peak_integral=ceil_h_step*sum(ext_peak(peak_points));

				temp_prof=ceil_bscatt(prof_num,bl_low_ind(nn):(c_h_ind-bl_high_ind(nn)));
				if max(temp_prof) <= bl_thres(nn)
					c_time_2(m)=c_time_1(n);
					h=c_heights_1(n);
					if nn == 1
						etaS=A*h^4+B*h^3+C*h^2+D*h+E;
					else
						etaS=h_etaS;
					end
					c_calfactor(m)=1/(2*etaS*peak_integral);
					m=m+1;
				end	
		
			end

			c_time_2=c_time_2(1:m-1);
			c_calfactor=c_calfactor(1:m-1);
			length(c_calfactor)
			cm=m-1;
			fprintf(fid_stats,'%d\t',cm);
			m=1;

			for kk=1:length(bin)
				indexbb=find(c_calfactor > bin(kk)-bin_step/2 & c_calfactor <= bin(kk)+bin_step/2);
				pdf(nn,kk)=pdf(nn,kk)+length(indexbb);
			end

			end	%time and height array lengths > 0
		end	%len1 > 0	
		nc=close(nc);



	end	%t(nn) > 0

end %for nn=1:4
max(max(pdf))
val.pdfout=pdf;


