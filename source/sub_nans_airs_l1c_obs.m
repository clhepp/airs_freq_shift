function [data_out] = sub_nans_airs_l1c_obs(data_in)

% function [] = sub_nans_airs_l1c_obs()
%
% Substitude NaNs in channels that are fake from AIRS L1b data
%
% Load the gap data:

% Check data_in
d1      = size(data_in);                    % data_in dimension list
[na nb] = size(data_in);
if(na ~= 2645) error('wrong input data, need 2645 channels'); whos data_in; return; end
data_in = reshape(data_in, na,nb);                       % make data_in a 2-D array


inFH = fopen('/home/chepplew/projects/airs/freqCal/data/gap4chan.asc','r');
  % skip 15 lines
  for i=1:15 junk = fgetl(inFH); end
  gap_array = textscan(inFH, '{%d {%d %d %d %d} %f %f %f}','Delimiter',',');
  l1c_gap_chans = gap_array{1};
fclose(inFH);

data_in(l1c_gap_chans,:) = NaN;
data_out = data_in;

% restore original shape of data_in
data_out = reshape(data_out, d1);

clear data_in;

%
