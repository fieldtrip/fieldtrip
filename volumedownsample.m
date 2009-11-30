function [down] = volumedownsample(cfg, source);

% VOLUMEDOWNSAMPLE downsamples an anatomical MRI or source reconstruction
% and optionally normalizes its coordinate axes, keeping the homogenous
% transformation matrix correct.
%
% Use as
%   [volume] = volumedownsample(cfg, volume)
% where the cconfiguration can contain
%   cfg.downsample = integer number (default = 1, i.e. no downsampling)
%   cfg.smooth     = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%
% This function is used by SOUREINTERPOLATE, SOURCEREAD and SOURCENORMALIZE.

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% check if SPM2 is in path and if not add
hastoolbox('SPM2',1);

%% checkdata see below!!! %%

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'unused',  {'voxelcoord'});

if ~isfield(cfg, 'downsample'), cfg.downsample = 1;     end
if ~isfield(cfg, 'keepinside'), cfg.keepinside = 'yes'; end
if ~isfield(cfg, 'parameter'),  cfg.parameter = 'all';  end
if ~isfield(cfg, 'smooth'),     cfg.smooth = 'no';      end

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter),
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

% check if the input data is valid for this function
source = checkdata(source, 'datatype', 'volume', 'feedback', 'no');

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

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
down.cfg = cfg;
