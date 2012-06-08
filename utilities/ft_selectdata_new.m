function [varargout] = ft_selectdata(cfg, varargin)

% FT_SELECTDATA makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% Use as
%  [data] = ft_selectdata_new(cfg, data, ...)

% Copyright (C) 2012, Robert Oostenveld
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

ft_defaults
ft_preamble help

% determine the characteristics of the input data
dtype = ft_datatype(varargin{1});
for i=2:length(varargin)
  % ensure that all subsequent inputs are of the same type
  ok = ft_datatype(varargin{i}, dtype);
  if ~ok, error('input data should be of the same datatype'); end
end

if strcmp(dtype, 'raw')
  
  % use selfromraw
  cfg.channel = ft_getopt(cfg, 'channel', 'all');
  cfg.latency = ft_getopt(cfg, 'latency', 'all');
  cfg.trials  = ft_getopt(cfg, 'trials',  'all');
  for i=1:length(varargin)
    varargin{i} = selfromraw(varargin{i}, 'rpt', cfg.trials, 'chan', cfg.channel, 'latency', cfg.latency);
  end
  
else
  
  hastime   = isfield(varargin{1}, 'time');
  hasfreq   = isfield(varargin{1}, 'freq');
  hasdimord = ~all(cellfun(@isempty, regexp(fieldnames(varargin{1}), '.*dimord')));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 2:
  %   ensure that the cfg is fully contained in the data and consistent over all inputs
  %   get the selection along each of the dimensions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if isfield(cfg, 'parameter') && isfield(varargin{1}, [cfg.parameter 'dimord'])
    dimord = [cfg.parameter 'dimord'];
  elseif isfield(varargin{1}, 'dimord')
    dimord = varargin{1}.dimord;
  end
  
  dimtok = tokenize(dimord, '_');
  dimsiz = nan(size(dimtok));
  dimfields = {};
  for i=1:numel(dimtok)
    % this switch-list is consistent with fixdimord
    switch dimtok{i}
      case 'time'
        dimsiz(i) = length(varargin{1}.time);
        dimfields{end+1} = 'time';
      case 'freq'
        dimsiz(i) = length(varargin{1}.freq);
        dimfields{end+1} = 'freq';
      case 'chan'
        dimsiz(i) = length(varargin{1}.label);
        dimfields{end+1} = 'label';
      case 'chancmb'
        dimsiz(i) = size(varargin{1}.labelcmb,1);
        dimfields{end+1} = 'labelcmb';
      case 'pos'
        dimsiz(i) = size(varargin{1}.pos,1);
        dimfields{end+1} = 'pos';
      case 'comp'
        dimsiz(i) = length(varargin{1}.label);
        dimfields{end+1} = 'label';
        
      case 'subj'
        error('FIXME');
      case 'rpt'
        error('FIXME');
      case 'rpttap'
        error('FIXME');
      case 'refchan'
        error('FIXME');
      case 'voxel'
        error('FIXME');
      case 'ori'
        error('FIXME');
      otherwise
        % try to guess the size from the corresponding field
        if isfield(varargin{1}, dimtok{i})
          siz = varargin{1}.(dimtok{i});
          if length(siz)==2 && any(siz==1)
            dimsiz(i) = prod(siz);
          end
        end
    end % switch
  end % for dimtok
  
  fn  = fieldnames(varargin{1})';
  sel = false(size(fn));
  for i=1:numel(fn)
    sel(i) = isequal(size(varargin{1}.(fn{i})), dimsiz);
  end
  
  % select the fields that represent the data
  datfields = fn(sel);
  
  switch dtype
    % this switch-list is consistent with ft_datatype
    
    case 'timelock'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
        [seltime, cfg] = getselection_time(cfg, varargin{i});
      end % varargin
      for i=1:numel(varargin)
        % get the selection from all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
        [seltime, cfg] = getselection_time(cfg, varargin{i});
        if ~isnan(selchan)
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan, datfields);
          varargin{i}.label = varargin{i}.label(selchan);
        end
        if ~isnan(seltime)
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime, datfields);
          varargin{i}.time  = varargin{i}.time(seltime);
        end
      end % varargin
      
    case 'freq'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
        [selfreq, cfg] = getselection_freq(cfg, varargin{i});
        if hastime
          [seltime, cfg] = getselection_time(cfg, varargin{i});
        end
      end % varargin
      keyboard
      
    case 'comp'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
      end % varargin
      keyboard
      
    case 'raw'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
      end % varargin
      keyboard
      
    case 'freqmvar'
      error('FIXME');
      
    case 'mvar'
      error('FIXME');
      
    case 'spike'
      error('FIXME');
      
    case 'volume'
      error('FIXME');
      
    case 'source'
      error('FIXME');
      
    case 'dip'
      error('FIXME');
      
    case 'chan'
      % this results from avgovertime/avgoverfreq after timelockstatistics or freqstatistics
      error('FIXME');
      
    otherwise
      % try to get the selection based on the field name
      seldim = cell(size(dimtok));
      for j=1:numel(seldim)
        seldim(j) = feval(['getselection_' dimtok{j}], cfg, varargin{i});
      end
  end
  
  % remove all fields from the data that do not pertain to the selection
  for i=1:numel(varargin)
    varargin{i} = keepfields(varargin{i}, [datfields dimfields {'cfg'}]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3:
%   if desired, concatenate over repetitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>2 && nargout==2
  % concatenate all input data into a single structure
  error('FIXME');
else
  % no reason to concatenate
  varargout = varargin;
  varargout{end+1} = cfg;
end

end % function ft_selectdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = keepfields(data, fn)

fn = setdiff(fieldnames(data), fn);
for i=1:numel(fn)
  data = rmfield(data, fn{i});
end

end % function keepfields

function [data] = makeselection(data, seldim, selindx, datfields)

  if numel(seldim) > 1
    for k = 1:numel(seldim)
      data = makeselection(data, seldim(k), selindx, datfields);
    end
  end

  for i=1:numel(datfields)
    switch seldim
      case 1
        data.(datfields{i}) = data.(datfields{i})(selindx,:,:,:,:,:);
      case 2
        data.(datfields{i}) = data.(datfields{i})(:,selindx,:,:,:,:);
      case 3
        data.(datfields{i}) = data.(datfields{i})(:,:,selindx,:,:,:);
      case 4
        data.(datfields{i}) = data.(datfields{i})(:,:,:,selindx,:,:);
      case 5
        data.(datfields{i}) = data.(datfields{i})(:,:,:,:,selindx,:);
      case 6
        data.(datfields{i}) = data.(datfields{i})(:,:,:,:,:,selindx);
      otherwise
        error('unsupported dimension (%d) for making a selection for %s', seldim, datfields{i});
    end % switch
  end % for datfields

end % function makeselection

function [chanindx, cfg] = getselection_chan(cfg, data)

% this return value specifies that no selection was specified
chanindx = nan;

if isfield(cfg, 'channel')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  chanindx = match_str(data.label, cfg.channel);
end

if isequal(chanindx, 1:length(data.label))
  % the cfg was updated, but no selection is needed for the data
  chanindx = nan;
end

end % function getselection_chan

function [timeindx, cfg] = getselection_time(cfg, data)
% possible specifications are
% cfg.latency = value     -> can be 'all'
% cfg.latency = [beg end] -> this is less common, preferred is to use toilim
% cfg.toilim  = [beg end]

assert(isfield(data, 'time'), 'the input data should have a time axis');

% this return value specifies that no selection was specified
timeindx = nan;

if isfield(cfg, 'latency')
  % deal with string selection
  if ischar(cfg.latency)
    if strcmp(cfg.latency, 'all')
      cfg.latency = [min(data.time) max(data.time)];
    else
      error('incorrect specification of cfg.latency');
    end
  end
  % deal with numeric selection
  if numel(cfg.latency)==1
    % this single value should be within the time axis of each input data structure
    tbin = nearest(data.time, cfg.latency, true, true);
    cfg.latency = data.time(tbin);
    timeindx = tbin;
  elseif numel(cfg.latency)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    mintime = min(data.time);
    maxtime = max(data.time);
    if all(cfg.latency<mintime) || all(cfg.latency>maxtime)
      error('the selected time range falls outside the time axis in the data');
    end
    tbeg = nearest(data.time, cfg.latency(1), false, false);
    tend = nearest(data.time, cfg.latency(2), false, false);
    cfg.latency = data.time([tbeg tend]);
    timeindx = tbeg:tend;
  elseif size(cfg.latency,2)==2
    % this may be used for specification of the computation, not for data selection
  else
    error('incorrect specification of cfg.latency');
  end
end % if cfg.latency

if isfield(cfg, 'toilim')
  if numel(cfg.toilim)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    mintime = min(data.time);
    maxtime = max(data.time);
    if all(cfg.toilim<mintime) || all(cfg.toilim>maxtime)
      error('the selected time range falls outside the time axis in the data');
    end
    tbin = nan(1,2);
    tbin(1) = nearest(data.time, cfg.toilim(1), false, false);
    tbin(2) = nearest(data.time, cfg.toilim(2), false, false);
    cfg.toilim = data.time(tbin);
  else
    error('incorrect specification of cfg.toilim');
  end
end % cfg.toilim

if isequal(timeindx, 1:length(data.time))
  % the cfg was updated, but no selection is needed for the data
  timeindx = nan;
end

end % function getselection_time

function [freqindx, cfg] = getselection_freq(cfg, data)
% possible specifications are
% cfg.frequency = value     -> can be 'all'
% cfg.frequency = [beg end] -> this is less common, preferred is to use foilim
% cfg.foilim    = [beg end]

assert(isfield(data, 'freq'), 'the input data should have a frequency axis');

% this return value specifies that no selection was specified
freqindx = nan;

if isfield(cfg, 'frequency')
  % deal with string selection
  if ischar(cfg.frequency)
    if strcmp(cfg.frequency, 'all')
      cfg.frequency = [min(data.freq) max(data.freq)];
    else
      error('incorrect specification of cfg.frequency');
    end
  end
  
  % deal with numeric selection
  for i=1:numel(varargin)
    if numel(cfg.latency)==1
      % this single value should be within the time axis of each input data structure
      fbin = nearest(data.freq, cfg.frequency, true, true);
      cfg.frequency = data.freq(fbin);
      freqindx = fbin;
    elseif numel(cfg.frequency)==2
      % the [min max] range can be specifed with +inf or -inf, but should
      % at least partially overlap with the freq axis of the input data
      minfreq = min(data.freq);
      maxfreq = max(data.freq);
      if all(cfg.frequency<minfreq) || all(cfg.frequency>maxfreq)
        error('the selected range falls outside the time axis in the data');
      end
      fbeg = nearest(data.freq, cfg.latency(1), false, false);
      fend = nearest(data.freq, cfg.latency(2), false, false);
      cfg.latency = data.freq([fbeg fend]);
      freqindx = fbeg:fend;
    elseif size(cfg.frequency,2)==2
      % this may be used for specification of the computation, not for data selection
    else
      error('incorrect specification of cfg.frequency');
    end
  end % for varargin
end % if cfg.frequency

if isfield(cfg, 'foilim')
  if numel(cfg.foilim)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    minfreq = min(data.freq);
    maxfreq = max(data.freq);
    if all(cfg.foilim<minfreq) || all(cfg.foilim>maxfreq)
      error('the selected range falls outside the frequency axis in the data');
    end
    fbin = nan(1,2);
    fbin(1) = nearest(data.freq, cfg.foilim(1), false, false);
    fbin(2) = nearest(data.freq, cfg.foilim(2), false, false);
    cfg.foilim = data.freq(fbin);
  else
    error('incorrect specification of cfg.foilim');
  end
end % cfg.foilim

if isequal(freqindx, 1:length(data.freq))
  % the cfg was updated, but no selection is needed for the data
  freqindx = nan;
end

end % function getselection_freq

function [rptindx, cfg] = getselection_rpt(cfg, data)
% this should deal with cfg.trials

% this return value specifies that no selection was specified
rptindx = nan;

error('FIXME');

end % function getselection_rpt

function [subjindx, cfg] = getselection_subj(cfg, data)
% this should deal with cfg.trials

% this return value specifies that no selection was specified
subjindx = nan;

error('FIXME');

end % function getselection_subj
