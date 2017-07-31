% check_airs_l1c_shift.m
%
% get original L1C granules

addpath /home/chepplew/myLib/matlib              % read_airs_l1c
addpath /asl/matlib/aslutil                      % rad2bt

warn_msgId = {'MATLAB:imagesci:hdf:removalWarningHDFSD',...
           'MATLAB:imagesci:hdfeos:removalWarningHDFSW'};
for i=1:length(warn_msgId) warning('off',warn_msgId{i});  end


% specify which channel to load (on the L1c grid)
achn  = 775; % (weak water line)
achn  = 284; % strong CO2 line
achns = [achn-9:achn+9];

% source directory
l1c_list = dir('/asl/data/airs/L1C/2013/002/AIRS.2013.01.02.*.hdf');
disp(['Found ' num2str(length(l1c_list)) ' granules']);

alat = []; alon = []; scntyp = []; ra = []; solzen = [];
for ifn = 1:length(l1c_list)
   %disp(['Loading file: ' l1c_list(ifn).name]);
   l1c      = read_airs_l1c([l1c_list(ifn).folder '/' l1c_list(ifn).name]);
   ra       = [ra; l1c.radiances(:,:,achns)];
   alat     = [alat; l1c.Latitude];
   alon     = [alon; l1c.Longitude];
   solzen   = [solzen; l1c.solzen];
   scntyp   = [scntyp, l1c.scan_node_type];
   granID   = unique(l1c.gindex);
   fprintf(1,'.')
end
fprintf(1,'\n');
dim1 = size(ra);
ra = permute(ra,[3,1,2]);
  whos alat alon solzen ra scntyp granID
clear l1c;

% Get drift corrected L1C granules

shift_list = dir('/home/chepplew/data/airs/L1C/2013/002/AIRS.2013.01.02.G*.mat');
disp(['Found ' num2str(length(shift_list)) ' shifted granules']);

ra2 = [];
for ifn = 1:length(shift_list)
  load([shift_list(ifn).folder '/' shift_list(ifn).name]);
  ra2   = [ra2; single(l1c.ra(:,:,achns))];
  fprintf(1,'.')
end
dim3  = size(ra2);
ra2   = permute(ra2,[3,1,2]);
dim2  = size(l1c.ra);
  whos ra ra2
  
% convert to BT
ba1 = real(rad2bt(l1c.freq(achns), ra));
ba2 = real(rad2bt(l1c.freq(achns), ra2));
bm1 = squeeze(nanmean(ba1,2));            bm1 = nanmean(bm1,2);
bm2 = squeeze(nanmean(ba2,2));            bm2 = nanmean(bm2,2);

figure(1);clf;plot(l1c.freq(achns),ba1(:,1,1),'.-', l1c.freq(achns),ba2(:,1,1),'.-');
figure(1);clf;plot(l1c.freq(achns),ba1(:,1,1) - ba2(:,1,1),'.-');
figure(1);clf;plot(l1c.freq(achns),bm1 - bm2,'.-');
% orbit variation of centre track FOV, 4th channel
ba1_cntr = ba1(:,:,45);
ba2_cntr = ba2(:,:,45);
xch = 12;   % 725.7247 cm-1.
figure(2);clf;plot(alat(:,45),ba1_cntr(xch,:),'.', alat(:,45),ba2_cntr(xch,:),'.');
figure(2);clf;plot(alat(:,45),ba1_cntr(xch,:) - ba2_cntr(xch,:),'.');
% subset day/night for FOV45
ind = find(solzen(:,45) < 90);      inn = find(solzen(:,45) >= 90);

figure(3);clf;
  h1=subplot(2,1,1);plot(alat(ind,45),ba1_cntr(xch,ind) - ba2_cntr(xch,ind),'.');
    grid on;ylabel('d(BT) K');title('2013.01.02 725.725wn shift');legend('day')
  h2=subplot(2,1,2);plot(alat(inn,45),ba1_cntr(xch,inn) - ba2_cntr(xch,inn),'.');
    grid on;ylabel('d(BT) K');xlabel('Latitude of Obs');legend('night')

  linkaxes([h1 h2],'x');set(h1,'xticklabel','');
  pp1=get(h1,'position');set(h1,'position',[pp1(1) pp1(2)-pp1(4)*0.1 pp1(3) pp1(4)*1.1])
  pp2=get(h2,'position');set(h2,'position',[pp2(1) pp2(2)+pp2(4)*0.1 pp2(3) pp2(4)*1.1])
  
