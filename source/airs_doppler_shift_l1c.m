function [shift_ppm, asc, l1c_list] = airs_doppler_shift_l1c(sdate)

% function [shift,asc,list] = airs_doppler_shift_l1c(date)
%
% SYNOPSIS: Returns doppler shifts for selected date from AIRS L1C granule data.
%
% INPUT: a datetime variable, for selected day to process. e.g. datetime(2013,01,02)
%
% OUTPUTS: shift: an array of doppler shifts for all 90 FOVs in ppm [90 x N] 
%          where N is the number of observations for the day - 2.
%          asc: a [1 x N] vector for ascending = 1, descending = 0.
%          l1c_list: strucutre array with fields referring to the L1C granules processed.
%
% NOTES:
% 1. the doppler shift calculation requires three adjacent satellite geolocation
%    values for fitting, therefore the first and last values for the day are sacrificed.
% 2.  Calls: calc_doppler_clh.m
% 3.  Dependencies: read_airs_l1c.m
%

warn_msg = 'MATLAB:imagesci:hdfeos:removalWarningHDFSW';
warning('off',warn_msg);

cd /home/chepplew/projects/airs/freqCal/run

addpath /home/chepplew/projects/airs/freqCal/source
addpath /home/chepplew/myLib/matlib                  % read_airs_l1c

% check input date
if(~isdatetime(sdate)) error('incorrect datetime value'); disp(sdate); return; end
disp(sdate);

% prepare date string for data file name
[nyr, nmon, nday] = datevec(sdate);
cyr = sprintf('%04d',nyr);  cmon = sprintf('%02d',nmon);  cday = sprintf('%02d',nday);
cdate = [cyr '.' cmon '.' cday];
jday  = datenum(sdate) - datenum(cyr,'YYYY') + 1;

l1c_list = dir(['/asl/data/airs/L1C/' cyr '/' sprintf('%03d',jday) ...
                '/AIRS.' cdate '*.hdf']);
disp(['Found ' num2str(numel(l1c_list)) ' L1C granule files']);

sat = struct; obs = struct;
sat_lat = []; sat_lon = []; sat_time = []; sat_hgt = []; 
obs_lat = []; obs_lon = []; obs_time = []; obs_topo = [];
node_type = [];  scan_node_type = [];
for i=1:length(l1c_list)
  l1c = read_airs_l1c([l1c_list(i).folder '/' l1c_list(i).name]);
  sat_lat    = [sat_lat,  l1c.sat_lat];
  sat_lon    = [sat_lon,  l1c.sat_lon];
  sat_time   = [sat_time, l1c.nadirTAI];
  sat_hgt    = [sat_hgt,  l1c.satheight];

  obs_lat    = [obs_lat;  l1c.Latitude];
  obs_lon    = [obs_lon;  l1c.Longitude];
  obs_time   = [obs_time; l1c.Time];
  obs_topo   = [obs_topo; l1c.topog];
  
  node_type(i,:) = l1c.node_type(1:10);
  scan_node_type = [scan_node_type; l1c.scan_node_type'];
  
  fprintf(1,'.');
end  
fprintf(1,'\n');
  whos obs_* sat_* node_type scan_node_type
nsat = size(sat_lat,2);              % also the number of cross-track scans.
[natrk nxtrk] = size(obs_lat);

% check natrk == nsat
if(natrk ~= nsat) error('unequal number of obs and satellite tracks'); return; end

% set up 3-point satellite vectors, and the Observation vectors
sat.lat1 = sat_lat(1:nsat-2);
sat.lat2 = sat_lat(2:nsat-1);
sat.lat3 = sat_lat(3:nsat);
sat.lon1 = sat_lon(1:nsat-2);
sat.lon2 = sat_lon(2:nsat-1);
sat.lon3 = sat_lon(3:nsat);
sat.hgt1 = 1000.*sat_hgt(1:nsat-2);
sat.hgt2 = 1000.*sat_hgt(2:nsat-1);
sat.hgt3 = 1000.*sat_hgt(3:nsat);
sat.tim1 = sat_time(1:nsat-2);
sat.tim2 = sat_time(2:nsat-1);
sat.tim3 = sat_time(3:nsat);
obs.lat  = obs_lat(2:nsat-1,:)';              % NB ' -> [1 x n]
obs.lon  = obs_lon(2:nsat-1,:)';
obs.top  = obs_topo(2:nsat-1,:)';
obs.tim  = obs_time(2:nsat-1,:)';
  obs 
  sat

shift_ppm = [];
for k = 1:90  
  obs.lat  = obs_lat(2:nsat-1,k)';              % NB ' -> [1 x n]
  obs.lon  = obs_lon(2:nsat-1,k)';
  obs.top  = obs_topo(2:nsat-1,k)';
  obs.tim  = obs_time(2:nsat-1,k)';
  % Call the Doppler calculator
  shift_ppm(k,:) = calc_doppler_clh(obs,sat);
end

asc = [];
for i = 2:numel(scan_node_type)-1             % to match computed doppler shift vectors
  scantyp = deblank(char(scan_node_type(i)));
  switch scantyp
    case 'A'
      asc(i-1) = 1;
    case 'D'
      asc(i-1) = 0;
  end
end


%{
% plot checks

figure(2);clf;plot(obs.lat,shift_ppm(1,:),'-', obs.lat,shift_ppm(45,:),'-', ...
   obs.lat,shift_ppm(90,:),'-');
inAsc = find(asc == 1);
inDsc = find(asc == 0);
figure(3);clf;plot(obs.lat(inAsc), shift_ppm(1,inAsc),'b.',...
  obs.lat(inDsc),shift_ppm(1,inDsc),'c.', ...
  obs.lat(inAsc),shift_ppm(90,inAsc),'r.',...
  obs.lat(inDsc),shift_ppm(90,inDsc),'m.');
  legend('FOV1.asc','FOV1.dsc','FOV90.asc','FOV90.dsc');
  xlabel('Observation Latitude'),ylabel('Doppler shift ppm');
  title('2013.01.02 AIRS doppler shift. All orbits');
  % saveas(gcf,'../figs/20130102_Doppler_FOV1_90_vs.lat.png','png')
  

   
%}
