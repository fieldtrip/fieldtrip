function [type, dimord] = ft_datatype(data, desired)

% FT_DATATYPE determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = ft_datatype(data)
%   [status]       = ft_datatype(data, desired)
%
% See also FT_DATATYPE_COMP FT_DATATYPE_FREQ FT_DATATYPE_MVAR
% FT_DATATYPE_SEGMENTATION FT_DATATYPE_PARCELLATION FT_DATATYPE_SOURCE
% FT_DATATYPE_TIMELOCK FT_DATATYPE_DIP FT_DATATYPE_HEADMODEL
% FT_DATATYPE_RAW FT_DATATYPE_SENS FT_DATATYPE_SPIKE FT_DATATYPE_VOLUME

% Copyright (C) 2008-2012, Robert Oostenveld
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

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip, segmentation, parcellation
israw          =  isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ~isfield(data,'trialtime');
isfreq         = (isfield(data, 'label') || isfield(data, 'labelcmb')) && isfield(data, 'freq') && ~isfield(data,'trialtime') && ~isfield(data,'origtrial'); %&& (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'cohspctrm') || isfield(data, 'fourierspctrm') || isfield(data, 'powcovspctrm'));
istimelock     =  isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ~isfield(data,'trialtime'); %&& ((isfield(data, 'avg') && isnumeric(data.avg)) || (isfield(data, 'trial') && isnumeric(data.trial) || (isfield(data, 'cov') && isnumeric(data.cov))));
iscomp         =  isfield(data, 'label') && isfield(data, 'topo') || isfield(data, 'topolabel');
isvolume       =  isfield(data, 'transform') && isfield(data, 'dim');
issource       =  isfield(data, 'pos');
isdip          =  isfield(data, 'dip');
ismvar         =  isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'lag'));
isfreqmvar     =  isfield(data, 'freq') && isfield(data, 'transfer');
ischan         =  isfield(data, 'dimord') && strcmp(data.dimord, 'chan') && ~isfield(data, 'time') && ~isfield(data, 'freq');
issegmentation = check_segmentation(data);
isparcellation = check_parcellation(data);

% check if it is a spike structure
spk_hastimestamp  = isfield(data,'label') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
spk_hastrials     = isfield(data,'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && isfield(data, 'trialtime') && isa(data.trialtime, 'numeric');
spk_hasorig       = isfield(data,'origtrial') && isfield(data,'origtime'); % for compatibility
isspike           = isfield(data, 'label') && (spk_hastimestamp || spk_hastrials || spk_hasorig);

if iscomp
  % comp should conditionally go before raw, otherwise the returned ft_datatype will be raw
  type = 'comp';
elseif israw
  type = 'raw';
elseif isfreqmvar
  % freqmvar should conditionally go before freq, otherwise the returned ft_datatype will be freq in the case of frequency mvar data
  type = 'freqmvar';
elseif isfreq
  type = 'freq';
elseif ismvar
  type = 'mvar';
elseif isdip
  % dip should conditionally go before timelock, otherwise the ft_datatype will be timelock
  type = 'dip';
elseif istimelock
  type = 'timelock';
elseif isspike
  type = 'spike';
elseif issegmentation
  % a segmentation data structure is a volume data structure, but in general not vice versa
  % segmentation should conditionally go before volume, otherwise the returned ft_datatype will be volume
  type = 'segmentation';
elseif isvolume
  type = 'volume';
elseif isparcellation
  % a parcellation data structure is a source data structure, but in general not vice versa
  % parcellation should conditionally go before source, otherwise the returned ft_datatype will be source
  type = 'parcellation';
elseif issource
  type = 'source';
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = check_segmentation(volume)
res = false;

if ~isfield(volume, 'dim')
  return
end

if isfield(volume, 'pos')
  return
end

if any(isfield(volume, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}))
  res = true;
  return
end

fn = fieldnames(volume);
for i=1:length(fn)
  if isfield(volume, [fn{i} 'label'])
    res = true;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = check_parcellation(source)
res = false;

if isfield(source, 'dim')
  return
end

if ~isfield(source, 'pos')
  return
end

if any(isfield(source, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}))
  res = true;
  return
end

fn = fieldnames(source);
for i=1:length(fn)
  if isfield(source, [fn{i} 'label'])
    res = true;
    return
  end
end
