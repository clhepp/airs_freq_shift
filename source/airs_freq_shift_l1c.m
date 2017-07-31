function [] = airs_freq_shift_l1c(sdate)

% AIRS frequency Correction in L1C data for single date.
%
% INPUT:  datetime variable for seelcted day to process.
%
% OUTPUT: None
%         Writes new L1C data file.
%
% Notes: to be used for AIRS L1c data
%        v1.0.0 is not stress tested and assumes whole granule files.
%
% Method:
%  1. Get the nominal L1b and L1c channel centers [2378] & [2645]
%  2. Sort channel sets, find common and synthetic channels in L1C set.
%  3. Load Observation data on L1c grid [2645] and compute orbit phase
%  4. Calculate frequency drift correction of channel centers [2378]
%  5. Resample radiances onto nominal grid and translate to L1C grid.
%  6. Save data to new L1C MAT files.
%
% v1.00 28Jul2017. C Hepplewhite. 

cd /home/chepplew/projects/airs/freqCal/run
addpath /home/chepplew/projects/airs/freqCal/source
addpath /home/chepplew/projects/airs/freqCal/data
addpath /home/chepplew/gitLib/JPL_nucal
addpath /home/chepplew/gitLib/l1c_nu_shift
addpath /home/chepplew/myLib/matlib                       % read_airs_l1c
addpath /asl/packages/ccast/source                        % seq_match.m
addpath /asl/matlib/h4tools
addpath /asl/packages/airs_decon/source

warn_msgId = {'MATLAB:imagesci:hdf:removalWarningHDFSD',...
           'MATLAB:imagesci:hdfeos:removalWarningHDFSW'};
for i=1:length(warn_msgId) warning('off',warn_msgId{i});  end

% Check input date
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

% ------------------------------------------------
% Directory to save data
% ------------------------------------------------
savD   = ['/home/chepplew/data/airs/L1C/' cyr '/' sprintf('%03d',jday) '/'];
if( exist(savD, 'dir') ~= 7) mkdir(savD); end

% load reference airs frequency grid:
load('airs_f.mat');
nL1b = length(f);                % 2378
nL1c = length(fairs);            % 2645

% Don't reload 'umbc_shift_l1c' parameters each time call function
persistent d1 
ix = ':';  
if isempty(d1)
    d1 = load('umbc_shift_1c');
    a = d1.a(ix); b = d1.b(ix);
end

% --------------------------------------------------
% Get cross-indexing of L1b and L1c channel sets
% --------------------------------------------------
% get the L1c channel set
%L1c_frq = fairs;      %load('freq2645.txt');
%nL1c    = length(L1c_frq);

% 1:2378 are the L1B channels, we want them sorted
[L1b_frq, ib] = sort(f);
[L1c_frq, ic] = sort(fairs);

% take the intersection of L1b and L1c channels [2314]
[ixb, ixc] = seq_match(L1b_frq, L1c_frq);

% any L1c channel not in the intersection is synthetic [331]
ind_syn = setdiff((1:nL1c)', ixc);
L1c_syn = L1c_frq(ind_syn);

% ------------------------------------------------
% Get v_nom on the L1b grid
% ------------------------------------------------
nom_date            = datetime(2010,01,22);
[l1b_nom_state]     = get_ab_state(nom_date);
yoff_nom            = get_yoff(nom_date);
[~,l1b_v_nom,~,~]   = gmodel(155.1325,yoff_nom,l1b_nom_state);

% put v_nom on the L1c grid
l1c_v_nom           = zeros(nL1c,1);
l1c_v_nom(ixc)      = l1b_v_nom(ixb);
l1c_v_nom(ind_syn)  = L1c_frq(ind_syn);

% --------------------------------------------------
% Get ab state and y_off for this date
% --------------------------------------------------
l1b_ab_sdate = get_ab_state(sdate);
yoff         = get_yoff(sdate);

% ---------------------------------------------------
% Get Observation data, One granule at a time
% ---------------------------------------------------
for ifn = 24:length(l1c_list)
  alat = []; alon = []; scntyp = []; ra = [];
  disp(['Loading file: ' l1c_list(ifn).name]);
  l1c      = read_airs_l1c([l1c_list(ifn).folder '/' l1c_list(ifn).name]);
  ra       = permute(l1c.radiances,[3 1 2]);
  alat     = [alat; l1c.Latitude];
  alon     = [alon; l1c.Longitude];
  scntyp   = [scntyp, l1c.scan_node_type];
  granID   = unique(l1c.gindex);

  dim_ra   = size(l1c.radiances);
    %whos alat alon ra scntyp granID
  
% Get orbit phase (0: desc eqX, 45: S.P., 90: asc eqX , 135: N.P., 180: desc eqX.)
% for interpolation record three values of phi (floor, actual and ceil)
  [na nx] = size(alat);
  phi  = [];
  for i = 1:na
    if( strcmp(char(scntyp(i)),'A') ) 
       phi(i,:) = [floor(0.5*(alat(i,45) + 180)), 0.5*(alat(i,45) + 180), ...
                   ceil(0.5*(alat(i,45) + 180))]; end
    if( strcmp(char(scntyp(i)),'D') ) 
       phi(i,:) = [floor(0.5*mod(360 - alat(i,45), 360)), 0.5*mod(360 - alat(i,45), 360), ...
                   ceil(0.5*(mod(360 - alat(i,45), 360)))]; end
  end

% --------------------------------------------
%       Channel drift correction
% -------------------------------------------
  disp('Calculating drift correction')
  yoff_a = []; yoff_b = []; yoff_p = []; new_freq = [];

% with interpolation - wrap at boundary 180+1 -> 1. (NB y=mod(x-1,180)+1 ensures y=[1:180])
  for i=1:na
    yoff_a = yoff(:,phi(i,1)+1); yoff_b = yoff(:,mod(phi(i,3),180)+1);
    for k = 1:17
      yoff_p(k,i) = interp1([phi(i,1)+1,phi(i,3)+1],[yoff_a(k),yoff_b(k)],phi(i,2)+1);
    end
    [f_lm,new_freq(:,i),m_lm,module] = gmodel(155.1325,yoff_p(:,i),l1b_ab_sdate);
    if(~mod(i,1000)) fprintf('.'); end
  end

%{
% without interpolation
  for i=1:na 
    yoff_p(:,i) = yoff(:,phi(i)+1);                         % NB +1 indexes from 1 to 180.
    [f_lm,new_freq(:,i),m_lm,module] = gmodel(155.1325,yoff_p(:,i),l1b_ab_sdate);
    if(~mod(i,1000)) fprintf('.'); end
  end 
%}  
% ----------------------------------------------
%   reconstitute L1c shift vector & resample
% ----------------------------------------------
  disp('Resampling radiances');
  ra_resamp = zeros(nL1c, na, nx);
  for i = 1:na
     l1c_shift           = zeros(nL1c,1);
     l1c_shift(ixc)      = new_freq(ixb,i);
     l1c_shift(ind_syn)  = L1c_frq(ind_syn);

   % Apply BT resample (on the L1C grid)
     abt     = real(rad2bt(l1c_v_nom, ra(:,i,:)));

     nu_in  = l1c_shift;
     dnu    = l1c_v_nom - nu_in;

     for j = 1:nx
       Tb_in     = abt(:,1,j);
       pp        = csape(nu_in, Tb_in);
       fprime    = fnder(pp);
       btderiv   = fnval(fprime,l1c_v_nom);
       Tb_resamp = Tb_in + (a.*btderiv + b).*dnu;

       ra_resamp(:,i,j) = bt2rad(l1c_v_nom,Tb_resamp);
     end
     fprintf(1,'.');
  end
  fprintf(1,'\n')

% granule file to save
  if(numel(granID) > 1) granID = granID(1); end
  savFn  = ['AIRS.' cdate '.' sprintf('G%03d',granID) '.L1C_shift.mat'];

% replace radiance and frequencies with new values
  l1c = rmfield(l1c,'radiances');
  l1c = rmfield(l1c,'freq');
  ra_resamp      = permute(ra_resamp, [2,3,1]);
  l1c.ra         = single(ra_resamp);
  l1c.freq       = l1c_v_nom;

% save to file:
  disp(['Saving file : ' savFn]);
  save(strcat(savD, savFn), 'l1c', '-v7.3');
  fprintf(1,'\n');
end                 % granule for-loop


% ---------------------- END OF FUNCTION -----------------------------
%{
% Plot checks
addpath /asl/matlib/aslutil
figure(1);clf;simplemap(alat(:,45), alon(:,45), phi);
figure(3);clf;plot([1:na], 1E6*(new_freq(640,:)-new_freq(640,1))./new_freq(640,1),'.-')
% substitude NaNs for synthetic channels
abt_nan       = sub_nans_airs_l1c_obs(abt);
ra_nan        = sub_nans_airs_l1c_obs(ra);
ra_resamp_nan = sub_nans_airs_l1c_obs(ra_resamp);

% check variation relative to nominal value of channel 640.
figure(2);clf;plot([1:size(freq,2)], freq(640,:)-l1b_v_nom(640),'.-')

% Check frequency spectrum of L1b+fake [2834]
figure(2);clf;plot([1:nL1f], sort(xfreq) - sort(xfreq_shift),'.');

% compare nominal freq spectrum as computed for 22Jan2010 and the f_airs:
figure(4);clf;plot([1:nL1c],fairs,'.', [1:nL1c],l1c_v_nom,'.');
figure(4);clf;plot([1:nL1c],fairs - l1c_v_nom,'.'); 

% check frequency spectrum of L1c [2645]
figure(2);clf;plot([1:nL1c],fairs,'.', [1:nL1c], l1c_shift,'.');
figure(2);clf;plot([1:nL1c], fairs - l1c_shift,'.');ylim([-0.05 0.05]);grid on;

% check a single spectrum shift values
figure(2);clf;plot([1:nL1c], 1E6*(fairs - l1c_shift)./fairs,'.');ylim([-15 15]);grid on;
 xlabel('L1C channel index');ylabel('relative shift ppm');title('Single spectrum shift values')

% check the resampled radiance vectors
abt        = real(rad2bt(l1c_v_nom, ra_nan));
abt_resamp = real(rad2bt(l1c_v_nom, ra_resamp_nan));

figure(4);clf;h1=subplot(2,1,1);plot(l1c_v_nom,abt(:,1,1),'-', l1c_v_nom,abt_resamp(:,1,1),'-'); 
  xlim([640 1150]);grid on;ylabel('BT K');
  h2=subplot(2,1,2);plot(l1c_v_nom,abt(:,1,1) - abt_resamp(:,1,1),'-'); 
  xlim([640 1150]);xlabel('wn cm-1');ylabel('nom - resamp K');


  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])

%}
