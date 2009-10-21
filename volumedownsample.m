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
% $Log: volumedownsample.m,v $
% Revision 1.26  2009/10/01 12:48:16  jansch
% allow for volumetric data with dimensionality > 3
%
% Revision 1.25  2009/05/19 15:22:30  jansch
% remove functional paramters prior to downsampling
%
% Revision 1.24  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.23  2009/01/12 13:05:20  sashae
% small change in call to checkconfig
%
% Revision 1.22  2008/11/21 13:56:12  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.21  2008/10/02 14:40:36  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.20  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.19  2008/09/17 14:53:34  roboos
% removed fixvolume (and underlying grid2transform), not needed any more because checkdata has the possibility of converting a pos to a transform
%
% Revision 1.18  2007/04/18 10:27:44  roboos
% no feedback for checkdata
%
% Revision 1.17  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.16  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.15  2006/09/13 09:49:15  roboos
% introduced(?) parameterselection function
%
% Revision 1.14  2006/08/16 10:52:53  marsie
% fixed use of spm_smooth
%
% Revision 1.13  2006/07/27 08:29:09  roboos
% updated documentation
%
% Revision 1.12  2006/07/27 08:01:53  roboos
% improved documentation
%
% Revision 1.11  2006/07/25 11:26:33  roboos
% fixed name of variable
%
% Revision 1.10  2006/07/19 12:22:53  jansch
% implemented smoothing for functional volumes
%
% Revision 1.9  2006/04/03 14:09:22  jansch
% removed conversion of inside-volume into vector of indices.
%
% Revision 1.8  2006/02/24 16:45:48  roboos
% switched to fixvolume function, with inside as 3d volume
% not neccessary any more to treat inside/outside seperately
% not neccessary any more to reshape volumes
%
% Revision 1.7  2006/01/31 09:57:38  jansch
% changed if-statement before reformatting the inside-volume into a vector,
% into explicitly checking for the presence of an inside-volume.
%
% Revision 1.6  2006/01/31 09:48:39  jansch
% fixed bug regarding last change
%
% Revision 1.5  2006/01/31 09:45:32  jansch
% changed else cfg.keepinside='no' into an elseif to prevent crashes
%
% Revision 1.4  2006/01/31 09:27:58  roboos
% in case of keepinside: instead of setting outside to [], remove it from structure
%
% Revision 1.3  2006/01/30 13:48:04  roboos
% switched order of parameterselection() and grid2transform() for consistency with other functions
%
% Revision 1.2  2006/01/24 14:20:35  roboos
% removed the obsolete option cfg.voxelcoord, new behaviour is that it is always 'yes'
%
% Revision 1.1  2006/01/05 12:58:19  roboos
% This function (VOLUMExxx) replaces a function with the name xxxVOLUME.
% The fields xgrid/ygrid/zgrid are removed (this is from now on handled by
% grid2transform and the VOLUMExxx function should only work on volumes that
% are described using a transformation matrix).
% Writing of spm/analyze volumes is handled by private/volumewrite_spm.
%
% Revision 1.9  2005/11/10 12:31:01  roboos
% prevent from downsampling if the requested downsampling is 1x
% this also prevents the double allocation of memory for the same data
%
% Revision 1.8  2005/09/29 00:35:20  roboos
% ensure that cfg.parameter is cell array
% add 'inside' to cfg.parameter if desired
%
% Revision 1.7  2005/08/19 16:05:08  roboos
% added default cfg.parameter=all
% switched to parameterselection() subfunction for looping over all interesting volumes
%
% Revision 1.6  2005/08/19 12:03:36  roboos
% add xgrid/ygrid/zgrid if not yet present
%
% Revision 1.5  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.4  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.3  2005/02/16 15:10:05  roboos
% added new volume fields from pcc beamformer
%
% Revision 1.2  2005/02/08 08:44:38  roboos
% fixed error in computation of new homogenous coordinate transformation matrix (in case of voxelcoords=yes), moved code to separate function (transformgrid)
% replaced computation for all scalar volumes by a for-loop over a parameter list
% implemented downsampling of inside/outside field in same way as the other volumes (keepinside=yes|no)
%
% Revision 1.1  2004/09/08 12:38:33  jansch
% Inserted version information, and moved function out of private folder. The
% function is the same to /private/downsamplevolume.m in older versions of the
% toolbox.
%
% Revision 1.2  2004/08/26 11:12:56  roboos
% added support for all known source parameters, updated help
%
% Revision 1.1  2004/08/26 09:12:26  roboos
% new implementation, to be used by different sourceXXX functions
%

fieldtripdefs

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
cfg.version.id = '$Id: volumedownsample.m,v 1.26 2009/10/01 12:48:16 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
down.cfg = cfg;
