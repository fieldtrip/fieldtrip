function comp = ft_datatype_comp(comp, varargin)

% FT_DATATYPE_COMP describes the FieldTrip MATLAB structure for comp data
%
% The comp data structure represents time-series channel-level data that has
% been decomposed or unmixed from the channel level into its components or
% "blind sources", for example using ICA (independent component analysis) or
% PCA. This data structure is usually generated with the FT_COMPONENTANALYSIS
% function.
%
% An example of a comp data structure with 149 components is shown here:
%
%          time: {1x10 cell}
%         trial: {1x10 cell}
%          topo: [149x149 double]
%     topolabel: {149x1 cell}
%         label: {149x1 cell}
%       fsample: 300
%           cfg: [1x1 struct]
%
% The only difference to the raw data structure is that the comp structure
% contains the additional fields topo and topolabel. See FT_DATATYPE_RAW for
% further details.
%
% Required fields:
%   - time, trial, label
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
% Revision history:
%
% (2003/latest) The initial version was defined
%
% See also FT_DATATYPE, FT_DATATYPE_RAW and FT_DATATYPE_xxx

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
version = keyval('version', varargin); if isempty(version), version = 'latest'; end

if strcmp(version, 'latest')
  version = '2003';
end

% convert it into a raw data structure
rawdata = comp;
rawdata = rmfield(rawdata, 'topo');
rawdata = rmfield(rawdata, 'topolabel');
rawdata = ft_datatype_raw(rawdata, 'version', version);

% add the component specific fields again
rawdata.topo      = comp.topo;
rawdata.topolabel = comp.topolabel;
comp = rawdata;
clear rawdata

