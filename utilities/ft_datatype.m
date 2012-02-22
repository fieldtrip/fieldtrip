function [type, dimord] = ft_datatype(data, desired)

% FT_DATATYPE determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = ft_datatype(data)
%   [type, dimord] = ft_datatype(data, desired)
%
% See also FT_CHANTYPE, FT_FILETYPE, FT_SENSTYPE, FT_VOLTYPE, FT_DATATYPE_COMP,
% FT_DATATYPE_DIP, FT_DATATYPE_FREQ, FT_DATATYPE_MVAR, FT_DATATYPE_RAW,
% FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE, FT_DATATYPE_TIMELOCK,
% FT_DATATYPE_VOLUME

% Copyright (C) 2008-2011, Robert Oostenveld
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

if isempty(data)
  return;
end

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      =  isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ~isfield(data,'trialtime');
isfreq     = (isfield(data, 'label') || isfield(data, 'labelcmb')) && isfield(data, 'freq') && ~isfield(data,'trialtime') && ~isfield(data,'origtrial'); %&& (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'cohspctrm') || isfield(data, 'fourierspctrm') || isfield(data, 'powcovspctrm'));
istimelock =  isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ~isfield(data,'trialtime'); %&& ((isfield(data, 'avg') && isnumeric(data.avg)) || (isfield(data, 'trial') && isnumeric(data.trial) || (isfield(data, 'cov') && isnumeric(data.cov))));
iscomp     =  isfield(data, 'label') && isfield(data, 'topo') || isfield(data, 'topolabel');
isvolume   =  isfield(data, 'transform') && isfield(data, 'dim');
issource   =  isfield(data, 'pos');
isdip      =  isfield(data, 'dip');
ismvar     =  isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'lag'));
isfreqmvar =  isfield(data, 'freq') && isfield(data, 'transfer');
ischan     =  isfield(data, 'dimord') && strcmp(data.dimord, 'chan') && ~isfield(data, 'time') && ~isfield(data, 'freq'); 
% check if isspike:
spk_hastimestamp = isfield(data,'label') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
spk_hastrials = isfield(data,'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ...
isfield(data, 'trialtime') && isa(data.trialtime, 'numeric');
spk_hasorig = isfield(data,'origtrial') && isfield(data,'origtime'); %% for compatibility
isspike = isfield(data, 'label') && (spk_hastimestamp || spk_hastrials || spk_hasorig);

if iscomp
  % comp should conditionally go before raw, otherwise the returned ft_datatype will be raw
  type = 'comp';  
elseif isfreqmvar
  % freqmvar should conditionally go before freq, otherwise the returned ft_datatype will be freq in the case of frequency mvar data
  type = 'freqmvar';
elseif ismvar
  type = 'mvar';
elseif israw
  type = 'raw';
elseif isfreq
  type = 'freq';
elseif istimelock
  type = 'timelock';
elseif isspike
  type = 'spike';
elseif isvolume
  type = 'volume';
elseif issource
  type = 'source';
elseif isdip
  type = 'dip';
elseif ischan
  % this results from avgovertime/avgoverfreq after timelockstatistics or freqstatistics
  type = 'chan';
else
  type = 'unknown';
end

if nargin>1
  % return a boolean value
  type = strcmp(type, desired);
  return;
end

if nargout>1
  % also return the dimord of the input data
  if isfield(data, 'dimord')
    dimord = data.dimord;
  else
    dimord = 'unknown';
  end
end

