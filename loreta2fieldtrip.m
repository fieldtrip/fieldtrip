function [source] = loreta2fieldtrip(filename, varargin)

% LORETA2FIELDTRIP reads and converts a LORETA source reconstruction into a
% FieldTrip data structure, which subsequently can be used for statistical
% analysis or other analysis methods implemented in Fieldtrip.
%
% Use as
%   [source]  =  loreta2fieldtrip(filename, ...)
% where optional arguments can be passed as key-value pairs.
%
% The following optional arguments are supported
%   'timeframe'  =  integer number, which timepoint to read (default is to read all)

% This function depends on the loreta_ind.mat file

% Copyright (C) 2006, Vladimir Litvak
%
% $Log: loreta2fieldtrip.m,v $
% Revision 1.3  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.2  2006/08/29 20:50:21  roboos
% major cleanup, added optional key-value inputs, implemented timeframe selection, implemented support for storing multiple timeframes in the source.avg.mom field (i.e. timecourses)
%
% Revision 1.1  2006/04/03 07:43:28  roboos
% first version under CVS control

fieldtripdefs

% get the optional input arguments
timeframe  =  keyval('timeframe', varargin); % will be empty if not specified

% start with an empty source structure
source  =  [];

if filetype(filename, 'loreta_slor')
  voxnumber    = 6239;
  lorind       = getfield(load('loreta_ind.mat'), 'ind_sloreta');
  source.dim   = size(lorind);
  source.xgrid =  -70:5:70;
  source.ygrid = -100:5:65;
  source.zgrid =  -45:5:70;
elseif filetype(filename, 'loreta_lorb')
  voxnumber    = 2394;
  lorind       = getfield(load('loreta_ind.mat'), 'ind_loreta');
  source.dim   = size(lorind);
  source.xgrid =  -66:7:67;
  source.ygrid = -102:7:66;
  source.zgrid =  -41:7:71;
else
  error('unsupported LORETA format');
end

source.transform = eye(4);      % FIXME the transformation matrix should be assigned properly
source.inside  = find(lorind ~= lorind(1));  % first voxel is outside
source.outside = find(lorind == lorind(1));  % first voxel is outside

fid = fopen(filename,'r', 'ieee-le');
% determine the length of the file
fseek(fid, 0, 'eof');
filesize = ftell(fid);
Ntime = filesize/voxnumber/4;

fprintf('file %s contains %d timepoints\n', filename, Ntime);
fprintf('file %s contains %d grey-matter voxels\n', filename, voxnumber);

if isempty(timeframe)
  % read the complete timecourses
  fseek(fid, 0, 'bof');
  activity = fread(fid, [voxnumber Ntime], 'float = >single');
elseif length(timeframe)==1
  % read only a single timeframe
  fseek(fid, 4*voxnumber*(timeframe-1), 'bof');
  activity = fread(fid, [voxnumber 1], 'float = >single');
else
  error('you can read either one timeframe, or the complete timecourse');
end

fclose(fid);

Ntime = size(activity,2);
if Ntime>1
  for i=1:voxnumber
    mom{i} = activity(i,:);
  end
  mom{end+1} = []; % this one is used
  source.avg.mom = mom(lorind);
  fprintf('returning the activity at %d timepoints as dipole moments for each voxel\n', Ntime);
else
  % put it in source.avg.pow
  activity(end+1) = nan;
  % reshuffle the activity to ensure that the ordering is correct
  source.avg.pow  = activity(lorind);
  fprintf('returning the activity at one timepoint as a single distribution of power\n');
end

% FIXME someone should figure out how to interpret the activity
fprintf('note that there is a discrepancy between dipole moment (amplitude) and power (amplitude squared)\n');

% add the options used here to the configuration
cfg = [];
cfg.timeframe = timeframe;
cfg.filename  = filename;
% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: loreta2fieldtrip.m,v 1.3 2008/09/22 20:17:43 roboos Exp $';
% remember the full configuration details
source.cfg = cfg;

