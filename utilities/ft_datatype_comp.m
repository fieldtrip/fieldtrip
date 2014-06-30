function comp = ft_datatype_comp(comp, varargin)

% FT_DATATYPE_COMP describes the FieldTrip MATLAB structure for comp data
%
% The comp data structure represents time-series channel-level data that has
% been decomposed or unmixed from the channel level into its components or
% "blind sources", for example using ICA (independent component analysis) or
% PCA. This data structure is usually generated with the FT_COMPONENTANALYSIS
% function.
%
% An example of a decomposed raw data structure with 100 components that resulted from
% a 151-channel MEG recording is shown here:
%
%       unmixing: [100x151 double]  the compoment unmixing matrix
%           topo: [151x100 double]  the compoment topographies
%      topolabel: {151x1 cell}      the channel labels (e.g. 'MRC13')
%           time: {1x10 cell}       the timeaxis [1*Ntime double] per trial
%          trial: {1x10 cell}       the numeric data [151*Ntime double] per trial
%          label: {100x1 cell}      the component labels (e.g. 'runica001')
%           grad: [1x1 struct]      information about the sensor array (for EEG it is called elec)
%            cfg: [1x1 struct]      the configuration used by the function that generated this data structure
%
% The only difference to the raw data structure is that the comp structure contains
% the additional fields unmixing, topo and topolabel. Besides representing the time
% series information as a raw data structure (see FT_DATATYPE_RAW), it is also
% possible for time series information to be represented as timelock or freq
% structures (see FT_DATATYPE_TIMELOCK or FT_DATATYPE_FREQ).
%
% Required fields:
%   - unmixing, topo, topolabel
%
% Optional fields:
%   - cfg, all fields from FT_DATATYPE_RAW, FT_DATATYPE_TIMELOCK or FT_DATATYPE_FREQ
%
% Historical fields:
%   - cfg, offset, fsample, grad, elec, label, sampleinfo, time, topo, topolabel, trial, unmixing, see bug2513
%
% Revision history:
% (2014) The combination of comp with raw, timelock or freq has been defined explicitly.
%
% (2011) The unmixing matrix has been added to the component data structure.
%
% (2003) The initial version was defined
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2011-2014, Robert Oostenveld
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
version       = ft_getopt(varargin, 'version', 'latest');
hassampleinfo = ft_getopt(varargin, 'hassampleinfo', []); % the default is determined in ft_datatype_raw
hastrialinfo  = ft_getopt(varargin, 'hastrialinfo', []);  % the default is determined in ft_datatype_raw

if strcmp(version, 'latest')
  version         = '2014';
  % the following are used further down
  rawversion      = 'latest';
  timelockversion = 'latest';
  freqversion     = 'latest';
end

if isempty(comp)
  return;
end

switch version
  case '2014'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(comp, 'unmixing')
      % in case the unmixing matrix is not present, construct the best estimate based on the mixing (topo) matrix
      % this is shared with the 2011 code
      if size(comp.topo,1)==size(comp.topo,2)
        comp.unmixing = inv(comp.topo);
      else
        comp.unmixing = pinv(comp.topo);
      end
    end
    
    % convert it into a raw data structure and update it to the latest version
    if ft_datatype(comp, 'raw')
      raw = comp;
      raw = rmfield(raw, 'topo');
      raw = rmfield(raw, 'unmixing');
      raw = rmfield(raw, 'topolabel');
      raw = ft_datatype_raw(raw, 'version', rawversion, 'hassampleinfo', hassampleinfo, 'hastrialinfo', hastrialinfo);
      
      % add the component specific fields again
      raw.unmixing  = comp.unmixing;
      raw.topo      = comp.topo;
      raw.topolabel = comp.topolabel;
      comp = raw;
      clear raw
      
    elseif ft_datatype(comp, 'timelock')
      timelock = comp;
      timelock = rmfield(timelock, 'topo');
      timelock = rmfield(timelock, 'unmixing');
      timelock = rmfield(timelock, 'topolabel');
      timelock = ft_datatype_timelock(timelock, 'version', timelockversion);
      
      % add the component specific fields again
      timelock.unmixing  = comp.unmixing;
      timelock.topo      = comp.topo;
      timelock.topolabel = comp.topolabel;
      comp = timelock;
      clear timelock
      
    elseif ft_datatype(comp, 'freq')
      freq = comp;
      freq = rmfield(freq, 'topo');
      freq = rmfield(freq, 'unmixing');
      freq = rmfield(freq, 'topolabel');
      freq = ft_datatype_freq(freq, 'version', freqversion);
      
      % add the component specific fields again
      freq.unmixing  = comp.unmixing;
      freq.topo      = comp.topo;
      freq.topolabel = comp.topolabel;
      comp = freq;
      clear freq
      
    end % raw, timelock or freq
    
  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(comp, 'unmixing')
      % in case the unmixing matrix is not present, construct the best estimate based on the mixing (topo) matrix
      % this is shared with the 2014 code
      if size(comp.topo,1)==size(comp.topo,2)
        comp.unmixing = inv(comp.topo);
      else
        comp.unmixing = pinv(comp.topo);
      end
    end
    
    if ft_datatype(comp, 'timelock') || ft_datatype(comp, 'freq')
      % timelock or freq were not supported in 2011, hence force conversion to raw data
      comp = ft_checkdata(comp, 'datatype', 'raw');
    end
    
  case '2003'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(comp, 'unmixing')
      % this field did not exist until November 2011
      comp = rmfield(comp, 'unmixing');
    end
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for comp datatype', version);
end


