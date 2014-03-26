function comp = ft_datatype_comp(comp, varargin)

% FT_DATATYPE_COMP describes the FieldTrip MATLAB structure for comp data
%
% The comp data structure represents time-series channel-level data that has
% been decomposed or unmixed from the channel level into its components or
% "blind sources", for example using ICA (independent component analysis) or
% PCA. This data structure is usually generated with the FT_COMPONENTANALYSIS
% function.
%
% An example of a comp data structure with 100 components that resulted from
% a 151-channel MEG recording is shown here:
%
%          time: {1x10 cell}
%         trial: {1x10 cell}
%      unmixing: [100x151 double]
%          topo: [151x100 double]
%     topolabel: {151x1 cell}
%         label: {100x1 cell}
%       fsample: 300
%           cfg: [1x1 struct]
%
% The only difference to the raw data structure is that the comp structure
% contains the additional fields unmixing, topo and topolabel. See
% FT_DATATYPE_RAW for further details.
%
% Required fields:
%   - time, trial, label, topo, unmixing
%
% Optional fields:
%   - sampleinfo, trialinfo, grad, elec, hdr, cfg
%
% Deprecated fields:
%   - fsample
%
% Obsoleted fields:
%   - offset
% 
% Historical fields:
%   - cfg, fsample, grad, label, sampleinfo, time, topo, topolabel, trial,
%   unmixing, see bug2513
%
% Revision history:
% (2011/latest) The unmixing matrix has been added to the component data
% structure.
%
% (2003) The initial version was defined
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
version       = ft_getopt(varargin, 'version', 'latest');
hassampleinfo = ft_getopt(varargin, 'hassampleinfo', []); % the default is determined in ft_datatype_raw
hastrialinfo  = ft_getopt(varargin, 'hastrialinfo', []);  % the default is determined in ft_datatype_raw

if strcmp(version, 'latest')
  compversion = '2011';
  rawversion  = 'latest';
else
  % Note that this does not ensure for backward compatibility support
  % that the exact old version of the comp structure will be recreated.
  % For example compversion=2007 will ensure fsample to be present, but
  % will not strip off the unmixing field
  rawversion = compversion;
end

if isempty(comp)
  return;
end

switch compversion
  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(comp, 'unmixing')
      % in case the unmixing matrix is not present, construct the best estimate
      % based on the mixing (topo) matrix
      if size(comp.topo,1)==size(comp.topo,2)
        comp.unmixing = inv(comp.topo);
      else
        comp.unmixing = pinv(comp.topo);
      end
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

% convert it into a raw data structure
rawdata = comp;
rawdata = rmfield(rawdata, 'topo');
rawdata = rmfield(rawdata, 'topolabel');
rawdata = ft_datatype_raw(rawdata, 'version', rawversion, 'hassampleinfo', hassampleinfo, 'hastrialinfo', hastrialinfo);

% add the component specific fields again
rawdata.unmixing  = comp.unmixing;
rawdata.topo      = comp.topo;
rawdata.topolabel = comp.topolabel;
comp = rawdata;
clear rawdata

