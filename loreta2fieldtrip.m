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
% 
% See also NUTMEG2FIELDTRIP, SPASS2FIELDTRIP, FIELDTRIP2SPSS,
% FT_SOURCEANALYSIS, FT_SOURCEPLOT

% This function depends on the loreta_ind.mat file

% Copyright (C) 2006, Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble callinfo

% get the optional input arguments
timeframe = ft_getopt(varargin, 'timeframe'); % will be empty if not specified

% start with an empty source structure
source  =  [];

if ft_filetype(filename, 'loreta_slor')
  voxnumber    = 6239;
  lorind       = getfield(load('loreta_ind.mat'), 'ind_sloreta');
  source.dim   = size(lorind);
  source.xgrid =  -70:5:70;
  source.ygrid = -100:5:65;
  source.zgrid =  -45:5:70;
elseif ft_filetype(filename, 'loreta_lorb')
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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble callinfo
ft_postamble history source

