function [down] = ft_volumedownsample(cfg, source);

% FT_VOLUMEDOWNSAMPLE downsamples an anatomical MRI or source reconstruction
% and optionally normalizes its coordinate axes, keeping the homogenous
% transformation matrix correct.
%
% Use as
%   [volume] = ft_volumedownsample(cfg, volume)
% where the cconfiguration can contain
%   cfg.downsample = integer number (default = 1, i.e. no downsampling)
%   cfg.smooth     = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%
% This function is used by FT_SOURCEINTERPOLATE, FT_VOLUMEWRITE and FT_VOLUMENORMALISE.

% Undocumented local options:
%   cfg.inputfile        = one can specifiy preanalysed saved data as input
%   cfg.outputfile       = one can specify output as file to save to disk

% Copyright (C) 2004, Robert Oostenveld
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

fieldtripdefs

%% ft_checkdata see below!!! %%

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'unused',  {'voxelcoord'});

if ~isfield(cfg, 'spmversion'), cfg.spmversion = 'spm8'; end
if ~isfield(cfg, 'downsample'), cfg.downsample = 1;     end
if ~isfield(cfg, 'keepinside'), cfg.keepinside = 'yes'; end
if ~isfield(cfg, 'parameter'),  cfg.parameter = 'all';  end
if ~isfield(cfg, 'smooth'),     cfg.smooth = 'no';      end
if ~isfield(cfg, 'inputfile'),  cfg.inputfile  = [];    end
if ~isfield(cfg, 'outputfile'), cfg.outputfile = [];    end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    source = loadvar(cfg.inputfile, 'source');
  end
end

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter),
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

% check if the input data is valid for this function
source = ft_checkdata(source, 'datatype', 'volume', 'feedback', 'no');

%make local copy of source and remove all functional parameters
param = parameterselection('all', source);
down  = source;
for k = 1:length(param)
  down = rmsubfield(down, param{k});
end

% select the parameters that should be downsampled
cfg.parameter = parameterselection(cfg.parameter, source);

% select the voxels that will be kept in the downsampled output volume
xsel = 1:cfg.downsample:source.dim(1);
ysel = 1:cfg.downsample:source.dim(2);
zsel = 1:cfg.downsample:source.dim(3);

% store the coordinate transformation and downsampled axes definition
down.transform = source.transform;
down.xgrid     = xsel;
down.ygrid     = ysel;
down.zgrid     = zsel;
down.dim = [length(xsel) length(ysel) length(zsel)];
if length(source.dim)>3,
  down.dim = [down.dim source.dim(4:end)];
end

% update the downsampled homogenous transformation matrix
down = grid2transform(down);

% smooth functional parameters, excluding anatomy and inside
if ~strcmp(cfg.smooth, 'no'),
  % check if the required spm is in your path:
  if strcmpi(cfg.spmversion, 'spm2'),
    ft_hastoolbox('SPM2',1);
  elseif strcmpi(cfg.spmversion, 'spm8'),
    ft_hastoolbox('SPM8',1);
  end

  for j = 1:length(cfg.parameter)
    if strcmp(cfg.parameter{j}, 'inside')
      fprintf('not smoothing %s\n', cfg.parameter{j});
    elseif strcmp(cfg.parameter{j}, 'anatomy')
      fprintf('not smoothing %s\n', cfg.parameter{j});
    else
      fprintf('smoothing %s with a kernel of %d voxels\n', cfg.parameter{j}, cfg.smooth);
      tmp = double(getsubfield(source, cfg.parameter{j}));
      spm_smooth(tmp, tmp, cfg.smooth);
      setsubfield(source, cfg.parameter{j}, tmp);
    end
  end
end

% downsample each of the parameters
if cfg.downsample~=1
  for i=1:length(cfg.parameter)
    fprintf('downsampling %s\n', cfg.parameter{i});
    tmp  = getsubfield(source, cfg.parameter{i});
    down = setsubfield(down, cfg.parameter{i}, tmp(xsel, ysel, zsel));    % downsample the volume
  end
else
  for i=1:length(cfg.parameter)
    fprintf('not downsampling %s\n', cfg.parameter{i});
    down = setsubfield(down, cfg.parameter{i}, reshape(getsubfield(source, cfg.parameter{i}), source.dim));
  end
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';
% remember the configuration details of the input data

try, cfg.previous = source.cfg; end

% remember the exact configuration details in the output
down.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'source', down); % use the variable name "data" in the output file
end

