function [downsample] = ft_volumedownsample(cfg, source)

% FT_VOLUMEDOWNSAMPLE downsamples an anatomical MRI or source reconstruction
% and optionally normalizes its coordinate axes, keeping the homogenous
% transformation matrix correct.
%
% Use as
%   [volume] = ft_volumedownsample(cfg, mri)
% where the input mri should be a single anatomical volume that was
% for example read with FT_READ_MRI or should be a volumetric source
% reconstruction resulting from FT_SOURCEANALYSIS or FT_SOURCEINTERPOLATE.
%
% The configuration can contain
%   cfg.downsample = integer number (default = 1, i.e. no downsampling)
%   cfg.parameter  = string, data field to downsample (default = 'all')
%   cfg.smooth     = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%   cfg.keepinside = 'yes' or 'no', keep the inside/outside labeling (default = 'yes')
%   cfg.spmversion = string, 'spm2' or 'spm8' (default = 'spm8')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEINTERPOLATE, FT_VOLUMEWRITE and FT_VOLUMENORMALISE

% Copyright (C) 2004-2014, Robert Oostenveld
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
ft_preamble init
ft_preamble debug
ft_preamble loadvar source
ft_preamble provenance source
ft_preamble trackconfig

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input data is valid for this function
source = ft_checkdata(source, 'datatype', 'volume', 'feedback', 'no');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'unused',  {'voxelcoord'});

if ~isfield(cfg, 'spmversion'), cfg.spmversion = 'spm8'; end
if ~isfield(cfg, 'downsample'), cfg.downsample = 1;      end
if ~isfield(cfg, 'keepinside'), cfg.keepinside = 'yes';  end
if ~isfield(cfg, 'parameter'),  cfg.parameter = 'all';   end
if ~isfield(cfg, 'smooth'),     cfg.smooth = 'no';       end

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter),
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

%make local copy of source and remove all functional parameters
param = parameterselection('all', source);
downsample = source;
for k = 1:length(param)
  downsample = rmsubfield(downsample, param{k});
end

% select the parameters that should be downsampled
cfg.parameter = parameterselection(cfg.parameter, source);

% select the voxels that will be kept in the downsampled output volume
xsel = 1:cfg.downsample:source.dim(1);
ysel = 1:cfg.downsample:source.dim(2);
zsel = 1:cfg.downsample:source.dim(3);

% store the coordinate transformation and downsampled axes definition
downsample.transform = source.transform;
downsample.xgrid     = xsel;
downsample.ygrid     = ysel;
downsample.zgrid     = zsel;
downsample.dim = [length(xsel) length(ysel) length(zsel)];
if length(source.dim)>3,
  downsample.dim = [downsample.dim source.dim(4:end)];
end

% update the downsampled homogenous transformation matrix
downsample = grid2transform(downsample);

% smooth functional parameters, excluding anatomy and inside
if isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no'),
  % check that SPM is on the path, try to add the preferred version
  if strcmpi(cfg.spmversion, 'spm2'),
    ft_hastoolbox('SPM2',1);
  elseif strcmpi(cfg.spmversion, 'spm8'),
    ft_hastoolbox('SPM8',1);
  elseif strcmpi(cfg.spmversion, 'spm12'),
    ft_hastoolbox('SPM12',1);
  end
  
  for j = 1:length(cfg.parameter)
    if strcmp(cfg.parameter{j}, 'inside')
      fprintf('not smoothing %s\n', cfg.parameter{j});
    elseif strcmp(cfg.parameter{j}, 'anatomy')
      fprintf('not smoothing %s\n', cfg.parameter{j});
    else
      tmp = volumesmooth(getsubfield(source, cfg.parameter{j}), cfg.smooth, cfg.parameter{j});
      setsubfield(source, cfg.parameter{j}, tmp);
    end
  end
end

% downsample each of the parameters
if cfg.downsample~=1
  for i=1:length(cfg.parameter)
    fprintf('downsampling %s\n', cfg.parameter{i});
    tmp        = getsubfield(source, cfg.parameter{i});
    downsample = setsubfield(downsample, cfg.parameter{i}, tmp(xsel, ysel, zsel));    % downsample the volume
  end
else
  for i=1:length(cfg.parameter)
    fprintf('not downsampling %s\n', cfg.parameter{i});
    downsample = setsubfield(downsample, cfg.parameter{i}, reshape(getsubfield(source, cfg.parameter{i}), source.dim));
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   source
ft_postamble provenance downsample
ft_postamble history    downsample
ft_postamble savevar    downsample
