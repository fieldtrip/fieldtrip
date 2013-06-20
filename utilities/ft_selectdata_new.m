function [varargout] = ft_selectdata_new(cfg, varargin)

% FT_SELECTDATA makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% Use as
%  [data] = ft_selectdata_new(cfg, data, ...)
%
% Valid cfg field are:
%   cfg.trials  = 1xN, trial indices to keep. It can be 'all'. You can use
%                 logical indexing, where false(1,N) removes all the trials) 
% For data with a time-dimension possible specifications are
%   cfg.latency = value     -> can be 'all'
%   cfg.latency = [beg end]
% For frequency data possible specifications are
%   cfg.frequency = value     -> can be 'all'
%   cfg.frequency = [beg end] -> this is less common, preferred is to use foilim
%   cfg.foilim    = [beg end]


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

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble help              % this will show the function help if nargin==0 and return an error
ft_preamble provenance        % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar varargin  % this reads the input data in case the user specified the cfg.inputfile option

% determine the characteristics of the input data
dtype = ft_datatype(varargin{1});
for i=2:length(varargin)
  % ensure that all subsequent inputs are of the same type
  ok = ft_datatype(varargin{i}, dtype);
  if ~ok, error('input data should be of the same datatype'); end
end

cfg = ft_checkconfig(cfg, 'renamed', {'toilim' 'latency'});

% this function only works for the new (2013x) source representation without sub-structures 
if strcmp(dtype, 'source')
  % update the old-style beamformer source reconstruction
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', '2013x');
  end
  if isfield(cfg, 'parameter') && length(cfg.parameter)>4 && strcmp(cfg.parameter(1:4), 'avg.')
    cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
  end
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
  
  cfg.trials  = ft_getopt(cfg, 'trials',  'all');
  if length(varargin)>1 && ~isequal(cfg.trials, 'all')
    error('it is ambiguous to a subselection of trials while concatenating data')
  end
  
  hastime   = isfield(varargin{1}, 'time');
  hasfreq   = isfield(varargin{1}, 'freq');
  hasdimord = ~all(cellfun(@isempty, regexp(fieldnames(varargin{1}), '.*dimord')));
  haspos    = isfield(varargin{1}, 'pos');
  
  avgoverchan = istrue(ft_getopt(cfg, 'avgoverchan', false));
  avgoverfreq = istrue(ft_getopt(cfg, 'avgoverfreq', false));
  avgovertime = istrue(ft_getopt(cfg, 'avgovertime', false));
  avgoverrpt  = istrue(ft_getopt(cfg, 'avgoverrpt',  false));
  
  % although being called region-of-interest, the selection is actually made over source positions
  avgoverpos  = istrue(ft_getopt(cfg, 'avgoverroi',  false));
  if avgoverpos
    for i=1:length(varargin)
      % must be a source representation, not a volume representation
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source');
    end
    % the rest of the implementation is not yet complete
    % there is already avgoverpos, which works, but the selection according to cfg.roi does not work
    error('this is not yet implemented, please see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2053')
  end
  
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
  dimfields = cell(size(dimtok));
  for i=1:numel(dimtok)
    % this switch-list is consistent with fixdimord
    switch dimtok{i}
      case 'time'
        dimsiz(i) = length(varargin{1}.time);
        dimfields{i} = 'time';
      case 'freq'
        dimsiz(i) = length(varargin{1}.freq);
        dimfields{i} = 'freq';
      case 'chan'
        dimsiz(i) = length(varargin{1}.label);
        dimfields{i} = 'label';
      case 'chancmb'
        dimsiz(i) = size(varargin{1}.labelcmb,1);
        dimfields{i} = 'labelcmb';
      case 'pos'
        dimsiz(i) = size(varargin{1}.pos,1);
        dimfields{i} = 'pos';
      case 'comp'
        dimsiz(i) = length(varargin{1}.label);
        dimfields{i} = 'label';
      case 'subj'
        % the number of elements along this dimension is implicit
        dimsiz(i) = nan;
        dimfields{i} = 'implicit';
      case 'rpt'
        % the number of elements along this dimension is implicit
        dimsiz(i) = nan;
        dimfields{i} = 'implicit';
      case 'rpttap'
        % the number of elements along this dimension is implicit
        dimsiz(i) = nan;
        dimfields{i} = 'implicit';
        
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
            dimfields{i} = dimtok{i};
          end
        end
    end % switch
  end % for dimtok
  dimfields{end+1} = 'dimord'; % also keep the dimord
  
  % deal with the data dimensions whose size is only implicitly represented
  if any(strcmp(dimfields, 'implicit'))
    fn  = fieldnames(varargin{1})';
    sel = false(size(fn));
    for i=1:numel(fn)
      if isequalwithoutnans(size(varargin{1}.(fn{i})), dimsiz)
        warning('using the "%s" field to determine the size along the unknown dimensions', fn{i});
        % update the size of all dimensions
        dimsiz = size(varargin{1}.(fn{i}));
        % update the fieldname of each dimension
        dimfields(strcmp(dimfields, 'implicit')) = dimtok(strcmp(dimfields, 'implicit'));
        break
      end
    end
    if any(strcmp(dimfields, 'implicit'))
      % it failed
      error('could not determine the size of the implicit "%s" dimension', dimfields{strcmp(dimfields, 'implicit')});
    end
  end
  
  fn  = fieldnames(varargin{1})';
  sel = false(size(fn));
  for i=1:numel(fn)
    sel(i) = isequal(size(varargin{1}.(fn{i})), dimsiz) || isequal(size(varargin{1}.(fn{i})), [dimsiz 1]);
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
        [selrpt,  cfg] = getselection_rpt (cfg, varargin{i}, 'datfields', datfields);
      end % varargin
      
      for i=1:numel(varargin)
        % get the selection from all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
        [seltime, cfg] = getselection_time(cfg, varargin{i});
        [selrpt,  cfg] = getselection_rpt (cfg, varargin{i}, 'datfields', datfields);
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan, avgoverchan, datfields);
        varargin{i} = makeselection_chan(varargin{i}, selchan, avgoverchan); % update the label field
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime, avgovertime, datfields);
        varargin{i} = makeselection_time(varargin{i}, seltime, avgovertime); % update the time field
        
        if ~any(isnan(selrpt))
          % FIXME could also be subject
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'rpt')), selrpt, avgoverrpt, datfields);
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
      
      for i=1:numel(varargin)
        % get the selection from all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
        [selfreq, cfg] = getselection_freq(cfg, varargin{i});
        [selrpt,  cfg, rptdim] = getselection_rpt (cfg, varargin{i}, 'datfields', datfields);
        if hastime
          [seltime, cfg] = getselection_time(cfg, varargin{i});
        end
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan, avgoverchan, datfields);
        varargin{i} = makeselection_chan(varargin{i}, selchan, avgoverchan); % update the label field
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'freq')), selfreq, avgoverfreq, datfields);
        varargin{i} = makeselection_freq(varargin{i}, selfreq, avgoverfreq); % update the freq field
        
        if ~any(isnan(selrpt))
          varargin{i} = makeselection(varargin{i}, rptdim, selrpt, avgoverrpt, datfields);
        end
        
        if hastime
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime, avgovertime, datfields);
          varargin{i} = makeselection_time(varargin{i}, seltime, avgovertime); % update the time field
        end
      end % varargin
      
    case 'comp'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
      end % varargin
      keyboard
      
    case 'raw'
      error('FIXME');
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selchan, cfg] = getselection_chan(cfg, varargin{i});
      end % varargin
      keyboard
      
    case 'source'
      for i=1:numel(varargin)
        % trim the selection to all inputs
        [selpos, cfg] = getselection_pos(cfg, varargin{i});
      end
      
      for i=1:numel(varargin)
        % get the selection from all inputs
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'pos')), selpos, avgoverpos, datfields);
        varargin{i} = makeselection_pos(varargin{i}, selpos, avgoverpos); % update the pos field
      end % varargin
      
    case 'freqmvar'
      error('FIXME');
      
    case 'mvar'
      error('FIXME');
      
    case 'spike'
      error('FIXME');
      
    case 'volume'
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
  end % switch dtype
  
  % update the dimord
  keep = {};
  if avgovertime
    dimtok = setdiff(dimtok, 'time');
  end
  if avgoverfreq
    dimtok = setdiff(dimtok, 'freq');
  end
  if avgoverpos
    dimtok = setdiff(dimtok, 'pos');
  else
    keep = [keep {'inside' 'outside' 'dim'}];
  end
  if avgoverrpt
    % FIXME could also be rpttap or subject
    dimtok = setdiff(dimtok, 'rpt');
  else
    keep = [keep {'cumtapcnt' 'cumsumcnt' 'sampleinfo' 'trialinfo'}];
  end
  for i=1:numel(varargin)
    varargin{i}.dimord = sprintf('%s_', dimtok{:});
    varargin{i}.dimord = varargin{i}.dimord(1:end-1);  % remove the last '_'
  end
  
  % remove all fields from the data that do not pertain to the selection
  for i=1:numel(varargin)
    varargin{i} = keepfields(varargin{i}, [datfields dimfields {'cfg' 'grad'} keep]);
  end
  
end % if raw or something else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3:
%   if desired, concatenate over repetitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout = varargin;

ft_postamble debug              % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig        % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance         % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous varargin  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history varargout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
% ft_postamble savevar varargout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option

if nargout>numel(varargout)
  % also return the input cfg
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

function data = makeselection(data, seldim, selindx, avgoverdim, datfields)

if numel(seldim) > 1
  for k = 1:numel(seldim)
    data = makeselection(data, seldim(k), selindx, datfields);
  end
end

for i=1:numel(datfields)
  if ~isnan(selindx)
    % the value NaN indicates that it is not needed to make a selection, rather take all values
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
  end
  if avgoverdim
    data.(datfields{i}) = mean(data.(datfields{i}), seldim);
  end
end % for datfields


end % function makeselection

function data = makeselection_chan(data, selchan, avgoverchan)
if avgoverchan && all(isnan(selchan))
  data.label = sprintf('%s+', data.label{:});        % concatenate all channel labels
  data.label = {data.label(1:end-1)};                  % remove the last '+'
elseif avgoverchan && ~any(isnan(selchan))
  data.label = sprintf('%s+', data.label{selchan});  % concatenate all channel labels
  data.label = {data.label(1:end-1)};                  % remove the last '+'
elseif ~isnan(selchan)
  data.label = data.label(selchan);
  data.label = data.label(:);
end
end % function makeselection_chan

function data = makeselection_freq(data, selfreq, avgoverfreq)
if avgoverfreq
  data = rmfield(data, 'freq');
elseif ~isnan(selfreq)
  data.freq  = data.freq(selfreq);
end
end % function makeselection_freq

function data = makeselection_time(data, seltime, avgovertime)
if avgovertime
  data = rmfield(data, 'time');
elseif ~isnan(seltime)
  data.time  = data.time(seltime);
end
end % function makeselection_time

function data = makeselection_pos(data, selpos, avgoverpos)
if avgoverpos
  data = rmfield(data, 'pos');
elseif ~isnan(selpos)
  data.pos = data.pos(selpos, :);
end
end % function makeselection_pos


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
% cfg.latency = [beg end]

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

% % Note: cfg.toilim handling removed as it was renamed to cfg.latency

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
  if numel(cfg.frequency)==1
    % this single value should be within the time axis of each input data structure
    if (data.freq ~= cfg.frequency)
      fbin = nearest(data.freq, cfg.frequency, true, true);
    else
      fbin = 1;
    end
    cfg.frequency = data.freq(fbin);
    freqindx = fbin;
  elseif numel(cfg.frequency)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the freq axis of the input data
    minfreq = min(data.freq);
    maxfreq = max(data.freq);
    if all(cfg.frequency<minfreq) || all(cfg.frequency>maxfreq)
      error('the selected range falls outside the frequency axis in the data');
    end
    fbeg = nearest(data.freq, cfg.frequency(1), false, false);
    fend = nearest(data.freq, cfg.frequency(2), false, false);
    cfg.frequency = data.freq([fbeg fend]);
    freqindx = fbeg:fend;
  elseif size(cfg.frequency,2)==2
    % this may be used for specification of the computation, not for data selection
  else
    error('incorrect specification of cfg.frequency');
  end
end % if cfg.frequency

if isfield(cfg, 'foilim')
  if ~ischar(cfg.foilim) && numel(cfg.foilim)==2
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

function [rptindx, cfg, rptdim] = getselection_rpt(cfg, data, varargin)
% this should deal with cfg.trials
datfields = ft_getopt(varargin, 'datfields');

if isfield(cfg, 'trials') && ~isequal(cfg.trials, 'all') && ~isempty(datfields)
  
  dimtok = tokenize(data.dimord, '_');
  rptdim = find(strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));
  
  if isempty(rptdim)
    % this return value specifies that no selection was specified
    rptindx = nan;
    return
  else
    rptindx = ft_getopt(cfg, 'trials');
    
    if islogical(rptindx)
      % convert from booleans to indices
      rptindx = find(rptindx);
    end
    
    rptindx = unique(sort(rptindx));
    rptindx = unique(sort(rptindx));
    rptsiz  = size(data.(datfields{1}), rptdim);
    
    if strcmp(dimtok{rptdim}, 'rpttap')
      %account for the tapers
      sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
      begtapcnt = sumtapcnt(1:end-1)+1;
      endtapcnt = sumtapcnt(2:end);
      begtapcnt = begtapcnt(rptindx);
      endtapcnt = endtapcnt(rptindx);
      tapers = zeros(1,sumtapcnt(end));
      for k = 1:length(begtapcnt)
        tapers(begtapcnt(k):endtapcnt(k)) = k;
      end
      rptindx   = find(tapers);
      [srt,ix] = sort(tapers(tapers~=0));
      rptindx  = rptindx(ix);
      %       cfg.trials = rptindx;
      % TODO FIXME think about whether this is a good or a bad thing...
      %warning('cfg.trials accounts for the number of tapers now');
    end
    
    if rptindx(1)<1
      error('cannot select rpt/subj/rpttap smaller than 1');
    elseif rptindx(end)>rptsiz
      error('cannot select rpt/subj/rpttap larger than the number of repetitions in the data');
    end
    
    % commented out because of rpttap dilemma...
    %     cfg.trials = rptindx;
    
    return
  end
  
else
  rptindx = nan;
  rptdim = nan;
end % if isfield cfg.trials

end % function getselection_rpt

function [posindx, cfg] = getselection_pos(cfg, data)
% possible specifications are <none>
posindx = 1:size(data.pos,1);
end % function getselection_pos

function ok = isequalwithoutnans(a, b)
% this is *only* used to compare matrix sizes, so we can ignore any
% singleton last dimension
numdiff = numel(b)-numel(a);

if numdiff > 0
  % assume singleton dimensions missing in a
  a = [a(:); ones(numdiff, 1)];
  b = b(:);
elseif numdiff < 0
  % assume singleton dimensions missing in b
  b = [b(:); ones(abs(numdiff), 1)];
  a = a(:);
end

c = ~isnan(a(:)) & ~isnan(b(:));
ok = isequal(a(c), b(c));

end % function isequalwithoutnans
