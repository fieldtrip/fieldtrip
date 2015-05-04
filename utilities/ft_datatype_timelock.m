function timelock = ft_datatype_timelock(timelock, varargin)

% FT_DATATYPE_TIMELOCK describes the FieldTrip MATLAB structure for timelock data
%
% The timelock data structure represents averaged or non-averaged
% event-releted potentials (ERPs, in case of EEG) or ERFs (in case
% of MEG). This data structure is usually generated with the
% FT_TIMELOCKANALYSIS or FT_TIMELOCKGRANDAVERAGE function.
%
% An example of a timelock structure containing the ERF for 151 channels
% MEG data is
%
%     dimord: 'chan_time'       defines how the numeric data should be interpreted
%        avg: [151x600 double]  the numeric data (in this example it contains the average values of the activity for 151 channels x 600 timepoints)  
%      label: {151x1 cell}      the channel labels (e.g. 'MRC13')
%       time: [1x600 double]    the timepoints in seconds
%        var: [151x600 double]  the variance of the activity for 151 channels x 600 timepoints
%       grad: [1x1 struct]      information about the sensor array (for EEG-data it is called elec)
%        cfg: [1x1 struct]      configuration structure used by the invoking FieldTrip function 
%
% Required fields:
%   - label, dimord, time
%
% Optional fields:
%   - avg, var, dof, cov, trial, grad, elec, cfg
%
% Deprecated fields:
%   - <none>
%
% Obsoleted fields:
%   - fsample
% 
% Historical fields:
%   - avg, cfg, cov, dimord, dof, dofvec, elec, fsample, grad, label,
% numcovsamples, numsamples, time, trial, var, see bug2513
%
% Revision history:
%
% (2011v2/latest) The description of the sensors has changed, see FT_DATATYPE_SENS
% for further information.
%
% (2011) The field 'fsample' was removed, as it was redundant.
%
% (2003) The initial version was defined.
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011, Robert Oostenveld
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

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2011v2';
end

if isempty(timelock)
  return;
end

% ensure consistency between the dimord string and the axes that describe the data dimensions
timelock = fixdimord(timelock);

% remove the unwanted fields, it is unclear when they were precisely used
if isfield(timelock, 'numsamples'),       timelock = rmfield(timelock, 'numsamples');       end
if isfield(timelock, 'numcovsamples'),    timelock = rmfield(timelock, 'numcovsamples');    end
if isfield(timelock, 'numblcovsamples'),  timelock = rmfield(timelock, 'numblcovsamples');  end

if ~iscolumn(timelock.label)
  timelock.label = timelock.label';
end
if ~isrow(timelock.time)
  timelock.time = timelock.time';
end
if ~isfield(timelock, 'label')
  warning('data structure is incorrect since it has no channel labels');
end

switch version
  case '2011v2'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(timelock, 'grad')
      % ensure that the gradiometer structure is up to date
      timelock.grad = ft_datatype_sens(timelock.grad);
    end
    
    if isfield(timelock, 'elec')
      % ensure that the electrode structure is up to date
      timelock.elec = ft_datatype_sens(timelock.elec);
    end

  case '2003'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are no known conversions for backward or forward compatibility support

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for timelock datatype', version);
end

