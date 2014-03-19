function [varargout] = ft_selectdata_new(cfg, varargin)

% FT_SELECTDATA makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% Use as
%  [data] = ft_selectdata_new(cfg, data, ...)
%
% The cfg artument is a configuration structure which can contain
%   cfg.channel = Nx1 cell-array with selection of channels (default = 'all'),
%                 see FT_CHANNELSELECTION
%   cfg.trials  = 1xN, trial indices to keep. It can be 'all'. You can use
%                 logical indexing, where false(1,N) removes all the trials) 
%
% For data with a time-dimension possible specifications are
%   cfg.latency = value     -> can be 'all'
%   cfg.latency = [beg end]
%
% For frequency data possible specifications are
%   cfg.frequency = value     -> can be 'all'
%   cfg.frequency = [beg end] -> this is less common, preferred is to use foilim
%   cfg.foilim    = [beg end]
%
% If multiple input arguments are provided, ft_selectdata will adjust
% the individual inputs such that either the intersection across inputs is
% retained (i.e. only the channel/time/frequency points that are shared
% across all input arguments), or the union across inputs is retained
% (replacing missing data with nans). In either case, the order (e.g. of
% the labels) is made consistent across inputs. Multiple inputs in combination
% with the selection of trials is not supported. The exact behavior can be
% specified with
%    cfg.selmode   = 'intersect' (default) or 'union'
%    cfg.tolerance = scalar, (default = 1e-5) tolerance value to determine 
%                    equality of time/frequency bins

% Copyright (C) 2012-2014, Robert Oostenveld & Jan-Mathijs Schoffelen
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
ft_preamble init              % this will reset warning_once and show the function help if nargin==0 and return an error
ft_preamble provenance        % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar varargin  % this reads the input data in case the user specified the cfg.inputfile option

cfg.tolerance = ft_getopt(cfg, 'tolerance', 1e-5);        % default tolerance for checking equality of time/freq axes
cfg.selmode   = ft_getopt(cfg, 'selmode',   'intersect'); % default is to take intersection, alternative 'union'

% determine the characteristics of the input data
dtype = ft_datatype(varargin{1});
for i=2:length(varargin)
  % ensure that all subsequent inputs are of the same type
  ok = ft_datatype(varargin{i}, dtype);
  if ~ok, error('input data should be of the same datatype'); end
end

cfg = ft_checkconfig(cfg, 'renamed', {'toilim' 'latency'});

% this function only works for the upcoming (not yet standard) source representation without sub-structures 
if strcmp(dtype, 'source')
  % update the old-style beamformer source reconstruction
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
  end
  if isfield(cfg, 'parameter') && length(cfg.parameter)>4 && strcmp(cfg.parameter(1:4), 'avg.')
    cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
  end
end

if strcmp(dtype, 'raw') || strcmp(dtype, 'comp')
  if strcmp(cfg.selmode, 'union')
    error('cfg.selmode ''union'' is not yet supported for raw data');
  end
    
  % use selfromraw
  cfg.channel = ft_getopt(cfg, 'channel', 'all', 1); % empty definition by user is meaningful
  cfg.latency = ft_getopt(cfg, 'latency', 'all', 1);
  cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1);
  
  for i=1:length(varargin)
    varargin{i} = selfromraw(varargin{i}, 'rpt', cfg.trials, 'chan', cfg.channel, 'latency', cfg.latency);
  end
  
else
  
  cfg.channel = ft_getopt(cfg, 'channel', 'all', 1);
  cfg.latency = ft_getopt(cfg, 'latency', 'all', 1);
  cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1);
  if ~isfield(cfg, 'foilim')
    cfg.frequency = ft_getopt(cfg, 'frequency', 'all', 1);
  end
  
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
  
  if strcmp(cfg.selmode, 'union') && (avgoverchan || avgoverfreq || avgovertime || avgoverrpt)
    error('cfg.selmode ''union'' in combination with averaging across one of the dimentions is not implemented');
  end
  
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
        fprintf('ft_selectdata: using the "%s" field to determine the size along the unknown dimensions\n', fn{i});
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
    sel(i) = (isequal(size(varargin{1}.(fn{i})), dimsiz)...
      || isequal(size(varargin{1}.(fn{i})), [dimsiz 1]))...
      && ~strcmp(fn{i}, 'label') && ~strcmp(fn{i}, 'time')...
      && ~strcmp(fn{i}, 'freq'); % make sure we do not treat a descriptive field as data
  end
  
  % select the fields that represent the data
  datfields = fn(sel);
  
  switch dtype
    % this switch-list is consistent with ft_datatype
    
    case 'timelock'
      % trim the selection to all inputs
      [selchan, cfg] = getselection_chan(cfg, varargin{:}, cfg.selmode);
      [seltime, cfg] = getselection_time(cfg, varargin{:}, cfg.tolerance, cfg.selmode);
      
      selrpt = cell(numel(varargin),1);
      for i=1:numel(varargin)
        [selrpt{i}] = getselection_rpt (cfg, varargin{i}, 'datfields', datfields);
         
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan{i}, avgoverchan, datfields, cfg.selmode);
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime{i}, avgovertime, datfields, cfg.selmode);
        % FIXME could also be subject
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'rpt')), selrpt{i}, avgoverrpt, datfields, 'intersect');
        
        varargin{i} = makeselection_chan(varargin{i}, selchan{i}, avgoverchan); % update the label field
        varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); % update the time field
        varargin{i} = makeselection_rpt (varargin{i}, selrpt{i}); % avgoverrpt for the supporting fields is dealt with later
        
        % make an exception for the covariance here (JM 20131128)
        if isfield(varargin{i}, 'cov') && (all(~isnan(selrpt{i})) || all(~isnan(selchan{i})))
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'chan'))+[0 1], selchan{i}, avgoverchan, {'cov'}, cfg.selmode);
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'rpt')),        selrpt{i},  avgoverrpt,  {'cov'}, 'intersect');
          datfields   = [datfields {'cov'}];
        end
        
        %shiftdim the datfields, because convention has it that this should
        %be done
        if avgoverrpt,
          for k =1:numel(datfields)
            varargin{i}.(datfields{k}) = shiftdim(varargin{i}.(datfields{k}),1);
          end
        end  
          
        
      end % varargin
      
      % in the case of selmode='union', create the union of the descriptive
      % axes
      if strcmp(cfg.selmode, 'union')
        label = varargin{1}.label;
        time = varargin{1}.time;
        
        for i=2:numel(varargin)
          tmplabel = varargin{i}.label;
          tmptime  = varargin{i}.time;
          
          time(~isfinite(time)) = tmptime(~isfinite(time));
         
          emptylabel = find(cellfun('isempty', label));
          for k=emptylabel(:)'
            label{k} = tmplabel{k};
          end
        end
        for i=1:numel(varargin)
          varargin{i}.label = label;
          varargin{i}.time  = time;
        end
      end
      
    case 'freq'
      % trim the selection to all inputs
      [selfreq, cfg] = getselection_freq(cfg, varargin{:}, cfg.tolerance, cfg.selmode);
      [selchan, cfg] = getselection_chan(cfg, varargin{:},                cfg.selmode); % tolerance not needed
      if hastime, [seltime, cfg] = getselection_time(cfg, varargin{:}, cfg.tolerance, cfg.selmode); end
      
      selrpt    = cell(numel(varargin),1);
      selrpttap = cell(numel(varargin),1);
      rptdim    = cell(numel(varargin),1);
      for i=1:numel(varargin)
        % the rpt selection stays within this loop, it only should work
        % with a single data argument anyway
        [selrpt{i}, ~, rptdim{i}, selrpttap{i}] = getselection_rpt(cfg, varargin{i}, 'datfields', datfields); % in case tapers were kept, selrpt~=selrpttap, otherwise selrpt==selrpttap
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan{i}, avgoverchan, datfields, cfg.selmode);
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'freq')), selfreq{i}, avgoverfreq, datfields, cfg.selmode);
        varargin{i} = makeselection(varargin{i}, rptdim{i}, selrpttap{i}, avgoverrpt, datfields, 'intersect');
        if hastime
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime{i}, avgovertime, datfields, cfg.selmode);
        end
      
        %shiftdim the datfields, because convention has it that this should
        %be done
        if avgoverrpt,
          for k =1:numel(datfields)
            varargin{i}.(datfields{k}) = shiftdim(varargin{i}.(datfields{k}),1);
          end
        end  
        
        varargin{i} = makeselection_chan(varargin{i}, selchan{i}, avgoverchan); % update the label field
        varargin{i} = makeselection_freq(varargin{i}, selfreq{i}, avgoverfreq); % update the freq field
        varargin{i} = makeselection_rpt(varargin{i}, selrpt{i}); % avgoverrpt for the supporting fields is dealt with later
        if hastime, 
          varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); % update the time field
        
          % also deal with the cumtapcnt-field, because it has a frequency
          % dimension when time dimension is present
          varargin{i} = makeselection_cumtapcnt(varargin{i}, selfreq{i}, avgoverfreq);
        end
      end % varargin
      
      % in the case of selmode='union', create the union of the descriptive
      % axes
      if strcmp(cfg.selmode, 'union')
        label = varargin{1}.label;
        freq  = varargin{1}.freq;
        if hastime,
          time = varargin{1}.time;
        end
        
        for i=2:numel(varargin)
          tmpfreq  = varargin{i}.freq;
          tmplabel = varargin{i}.label;
          if hastime,
            tmptime  = varargin{i}.time;
          end
          
          freq(~isfinite(freq)) = tmpfreq(~isfinite(freq));
          if hastime,
            time(~isfinite(time)) = tmptime(~isfinite(time));
          end
          
          emptylabel = find(cellfun('isempty', label));
          for k=emptylabel(:)'
            label{k} = tmplabel{k};
          end
        end
        for i=1:numel(varargin)
          varargin{i}.freq  = freq;
          varargin{i}.label = label;
          if hastime,
            varargin{i}.time = time;
          end
        end
      end
      
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
    [dum, order] = setdiff(dimtok, 'time');
    dimtok       = dimtok(order);
  end
  if avgoverfreq
    [dum, order] = setdiff(dimtok, 'freq');
    dimtok       = dimtok(order);
  end
  if avgoverpos
    [dum, order] = setdiff(dimtok, 'pos');
    dimtok       = dimtok(order);  
  else
    keep = [keep {'inside' 'outside' 'dim'}];
  end
  if avgoverrpt
    sel = ismember(dimtok, {'rpt' 'rpttap' 'subject'});
    [dum, order] = setdiff(dimtok, dimtok{sel});
    dimtok       = dimtok(order);
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
% ft_postamble previous varargin  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
% ft_postamble history varargout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg

% note that the cfg.previous thingy does not work with the postamble,
% because the postamble puts the cfgs of all input arguments in the (first)
% output argument's xxx.cfg
for k = 1:numel(varargout)
  varargout{k}.cfg          = cfg;
  if isfield(varargin{k}, 'cfg')
    varargout{k}.cfg.previous = varargin{k}.cfg;
  end
end

% ft_postamble savevar varargout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option

if nargout>numel(varargout)
  % also return the input cfg
  varargout{end+1} = cfg;
end

end % function ft_selectdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = keepfields(data, fn)

% KEEPFIELDS behaves opposite to RMFIELD: it keeps the specified
% fieldnames, and removes the rest.
%
% Use as
%   data = keepfields(data, fn)
%
% Where fn is a cell-array of fields to keep

fn = setdiff(fieldnames(data), fn);
for i=1:numel(fn)
  data = rmfield(data, fn{i});
end

end % function keepfields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection(data, seldim, selindx, avgoverdim, datfields, selmode)

if numel(seldim) > 1
  for k = 1:numel(seldim)
    data = makeselection(data, seldim(k), selindx, avgoverdim, datfields, selmode);
  end
  return;
end

switch selmode
  case 'intersect'
    
    for i=1:numel(datfields)
      
      if isempty(selindx) || all(~isnan(selindx))
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
    
  case 'union'
    
    for i=1:numel(datfields)
      tmp = data.(datfields{i});
      siz = size(tmp);
      siz(seldim) = numel(selindx);
      data.(datfields{i}) = nan+zeros(siz);
      sel = isfinite(selindx);
      switch seldim
        case 1
          data.(datfields{i})(sel,:,:,:,:,:) = tmp(selindx(sel),:,:,:,:,:);
        case 2
          data.(datfields{i})(:,sel,:,:,:,:) = tmp(:,selindx(sel),:,:,:,:);
        case 3
          data.(datfields{i})(:,:,sel,:,:,:) = tmp(:,:,selindx(sel),:,:,:);
        case 4
          data.(datfields{i})(:,:,:,sel,:,:) = tmp(:,:,:,selindx(sel),:,:);
        case 5
          data.(datfields{i})(:,:,:,:,sel,:) = tmp(:,:,:,:,selindx(sel),:);
        case 6
          data.(datfields{i})(:,:,:,:,:,sel) = tmp(:,:,:,:,:,selindx(sel));
        otherwise
          error('unsupported dimension (%d) for making a selection for %s', seldim, datfields{i});
      end
    end
    if avgoverdim
      data.(datfields{i}) = mean(data.(datfields{i}), seldim);
    end
    
  otherwise
end

end % function makeselection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_chan(data, selchan, avgoverchan)
if avgoverchan && all(isnan(selchan))
  str = sprintf('%s, ', data.label{:});
  str = str(1:(end-2));
  str = sprintf('mean(%s)', str);
  data.label = {str};
elseif avgoverchan && ~any(isnan(selchan))
  str = sprintf('%s, ', data.label{selchan});
  str = str(1:(end-2));
  str = sprintf('mean(%s)', str);
  data.label = {str};                 % remove the last '+'
elseif all(isfinite(selchan))
  data.label = data.label(selchan);
  data.label = data.label(:);
elseif numel(selchan)==1 && any(~isfinite(selchan))
  % do nothing
elseif numel(selchan)>1  && any(~isfinite(selchan))
  tmp = cell(numel(selchan),1);
  for k = 1:numel(tmp)
    if isfinite(selchan(k))
      tmp{k} = data.label{selchan(k)};
    end
  end
  data.label = tmp;
elseif isempty(selchan)
  data.label = {};
end
end % function makeselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_freq(data, selfreq, avgoverfreq)
if avgoverfreq
  %data = rmfield(data, 'freq');
  if ~isnan(selfreq)
    data.freq  = mean(data.freq(selfreq));
  else
    data.freq  = mean(data.freq);
  end
elseif numel(selfreq)==1 && ~isfinite(selfreq)
  % do nothing
elseif numel(selfreq)==1 && isfinite(selfreq)
  data.freq = data.freq(selfreq);
elseif numel(selfreq)>1 && any(~isfinite(selfreq))
  tmp = selfreq(:)';
  sel = isfinite(selfreq);
  tmp(sel)  = data.freq(selfreq(sel));
  data.freq = tmp;
elseif numel(selfreq)>1 && all(isfinite(selfreq))
  data.freq = data.freq(selfreq);
elseif isempty(selfreq)
  data.freq  = zeros(1,0);
end
end % function makeselection_freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_cumtapcnt(data, selfreq, avgoverfreq)

if ~isfield(data, 'time')
  error('the subfunction makeselection_cumtapcnt should only be called when there is a time dimension in the data');
end
if ~isfield(data, 'cumtapcnt')
  return;
end

if avgoverfreq
  %data = rmfield(data, 'freq');
  if ~isnan(selfreq)
    data.cumtapcnt  = mean(data.cumtapcnt(:,selfreq),2);
  else
    data.cumtapcnt  = mean(data.cumtapcnt,2);
  end
elseif numel(selfreq)==1 && ~isfinite(selfreq)
  % do nothing
elseif numel(selfreq)==1 && isfinite(selfreq)
  data.cumtapcnt = data.cumtapcnt(:,selfreq);
elseif numel(selfreq)>1 && any(~isfinite(selfreq))
  tmp  = selfreq(:)';
  tmp2 = zeros(size(data.cumtapcnt,1), numel(selfreq));
  sel = isfinite(selfreq);
  tmp2(:, sel)  = data.cumtapcnt(:,selfreq(sel));
  data.freq = tmp2;
elseif numel(selfreq)>1 && all(isfinite(selfreq))
  data.cumtapcnt = data.cumtapcnt(:,selfreq);
elseif isempty(selfreq)
  %data.cumtapcnt  = zeros(1,0);
end
end % function makeselection_cumtapcnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_time(data, seltime, avgovertime)
if avgovertime
  data = rmfield(data, 'time');
elseif numel(seltime)==1 && ~isfinite(seltime)
  % do nothing
elseif numel(seltime)==1 && isfinite(seltime)
  data.time = data.time(seltime);
elseif numel(seltime)>1 && any(~isfinite(seltime))
  tmp = seltime(:)';
  sel = isfinite(seltime);
  tmp(sel)  = data.time(seltime(sel));
  data.time = tmp;
elseif numel(seltime)>1 && all(isfinite(seltime))
  data.time = data.time(seltime);
elseif isempty(seltime)
  data.time  = zeros(1,0);
end
end % function makeselection_time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_pos(data, selpos, avgoverpos)
if avgoverpos
  data = rmfield(data, 'pos');
elseif ~isnan(selpos)
  data.pos = data.pos(selpos, :);
end
end % function makeselection_pos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_rpt(data, selrpt)
if all(isfinite(selrpt)) || isempty(selrpt)
  if isfield(data, 'cumtapcnt')
    data.cumtapcnt = data.cumtapcnt(selrpt,:,:);
  end
  if isfield(data, 'cumsumcnt')
    data.cumsumcnt = data.cumsumcnt(selrpt,:,:);
  end
  if isfield(data, 'trialinfo')
    data.trialinfo = data.trialinfo(selrpt,:);
  end
  if isfield(data, 'sampleinfo')
    data.sampleinfo = data.sampleinfo(selrpt,:);
  end
end
end % function makeselection_rpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chanindx, cfg] = getselection_chan(cfg, varargin)

ndata = numel(varargin)-1;
selmode = varargin{end};

% loop over data once to initialize 
chanindx = cell(numel(varargin)-1,1);
label    = cell(1,0);
if isfield(cfg, 'channel')

  for k = 1:ndata
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label);
    label      = union(label, selchannel);
  end
  
  indx = nan+zeros(numel(label), ndata);

  for k = 1:ndata
    [ix, iy] = match_str(label, varargin{k}.label);
    indx(ix,k) = iy;
  end

  switch selmode
    case 'intersect'
      sel      = sum(isfinite(indx),2)==ndata;
      indx     = indx(sel,:);
      label    = varargin{1}.label(indx(:,1));
    case 'union'
      % don't do a subselection
    otherwise
  end

  for k = 1:ndata
    chanindx{k,1} = indx(:,k);
  end
  cfg.channel = label;
else
  
  for k = 1:ndata
    chanindx{k,1} = nan;
  end
  
end

end % function getselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeindx, cfg] = getselection_time(cfg, varargin)

% possible specifications are
% cfg.latency = value     -> can be 'all'
% cfg.latency = [beg end]
ndata = numel(varargin)-2;
tol   = varargin{end-1};
selmode = varargin{end};

% loop over data once to initialize 
timeindx = cell(numel(varargin)-2,1);
timeaxis = zeros(1,0);
for k = 1:ndata
  assert(isfield(varargin{k}, 'time'), 'the input data should have a time axis');

  % this return value specifies that no selection was specified
  timeindx{k,1} = nan;
  
  % update the axis along which the frequencies are defined
  timeaxis = union(timeaxis, round(varargin{k}.time(:)/tol)*tol); 
end

indx = nan+zeros(numel(timeaxis), ndata);
for k = 1:ndata
  [~, ix, iy] = intersect(timeaxis, round(varargin{k}.time(:)/tol)*tol);
  indx(ix,k) = iy;
end

switch selmode
  case 'intersect'
    sel      = sum(isfinite(indx),2)==ndata;
    indx     = indx(sel,:);
    timeaxis = varargin{1}.time(indx(:,1));
  case 'union'
    % don't do a subselection
  otherwise
end

if isfield(cfg, 'latency')
  % deal with string selection
  if ischar(cfg.latency)
    if strcmp(cfg.latency, 'all')
      cfg.latency = [min(timeaxis) max(timeaxis)];
    else
      error('incorrect specification of cfg.latency');
    end
  end
  % deal with numeric selection
  if numel(cfg.latency)==1
    % this single value should be within the time axis of each input data structure
    tbin = nearest(timeaxis, cfg.latency, true, true);
    cfg.latency = timeaxis(tbin);
    
    for k = 1:ndata
      timeindx{k,1} = indx(tbin, k);
    end
    
  elseif numel(cfg.latency)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    mintime = min(timeaxis);
    maxtime = max(timeaxis);
    if all(cfg.latency<mintime) || all(cfg.latency>maxtime)
      error('the selected time range falls outside the time axis in the data');
    end
    tbeg = nearest(timeaxis, cfg.latency(1), false, false);
    tend = nearest(timeaxis, cfg.latency(2), false, false);
    cfg.latency = timeaxis([tbeg tend]);
    
    for k = 1:ndata
      timeindx{k,1} = indx(tbeg:tend, k);
    end
    
  elseif size(cfg.latency,2)==2
    % this may be used for specification of the computation, not for data selection
  elseif isempty(cfg.latency)
    
    for k = 1:ndata
      timeindx{k,1} = [];
    end
  else
    error('incorrect specification of cfg.latency');
  end
end % if cfg.latency

% % Note: cfg.toilim handling removed as it was renamed to cfg.latency
for k = 1:ndata
  if isequal(timeindx, 1:length(timeaxis))
    % the cfg was updated, but no selection is needed for the data
    timeindx{k,1} = nan;
  end
end

end % function getselection_time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freqindx, cfg] = getselection_freq(cfg, varargin)
% possible specifications are
% cfg.frequency = value     -> can be 'all'
% cfg.frequency = [beg end] -> this is less common, preferred is to use foilim
% cfg.foilim    = [beg end]

ndata = numel(varargin)-2;
tol   = varargin{end-1};
selmode = varargin{end};

% loop over data once to initialize 
freqindx = cell(numel(varargin)-2,1);
freqaxis = zeros(1,0);
for k = 1:ndata
  assert(isfield(varargin{k}, 'freq'), 'the input data should have a frequency axis');

  % this return value specifies that no selection was specified
  freqindx{k,1} = nan;
  
  % update the axis along which the frequencies are defined
  freqaxis = union(freqaxis, round(varargin{k}.freq(:)/tol)*tol); 
end

indx = nan+zeros(numel(freqaxis), ndata);
for k = 1:ndata
  [~, ix, iy] = intersect(freqaxis, round(varargin{k}.freq(:)/tol)*tol);
  indx(ix,k) = iy;
end

switch selmode
  case 'intersect'
    sel      = sum(isfinite(indx),2)==ndata;
    indx     = indx(sel,:);
    freqaxis = varargin{1}.freq(indx(:,1));
  case 'union'
    % don't do a subselection
  otherwise
end

if isfield(cfg, 'frequency')
  % deal with string selection
  if ischar(cfg.frequency)
    if strcmp(cfg.frequency, 'all')
      cfg.frequency = [min(freqaxis) max(freqaxis)];
    else
      error('incorrect specification of cfg.frequency');
    end
  end
  
  % deal with numeric selection
  if numel(cfg.frequency)==1
    % this single value should be within the frequency axis of each input data structure
    fbin = nearest(freqaxis, cfg.frequency, true, true);
    cfg.frequency = freqaxis(fbin);
    
    for k = 1:ndata
      freqindx{k,1} = indx(fbin,k);
    end
    
  elseif numel(cfg.frequency)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the freq axis of the input data
    minfreq = min(freqaxis);
    maxfreq = max(freqaxis);
    if all(cfg.frequency<minfreq) || all(cfg.frequency>maxfreq)
      error('the selected range falls outside the frequency axis in the data');
    end
    fbeg = nearest(freqaxis, cfg.frequency(1), false, false);
    fend = nearest(freqaxis, cfg.frequency(2), false, false);
    cfg.frequency = freqaxis([fbeg fend]);
    
    for k = 1:ndata
      freqindx{k,1} = indx(fbeg:fend,k);
    end
    
  elseif size(cfg.frequency,2)==2
    % this may be used for specification of the computation, not for data selection
  elseif isempty(cfg.frequency)
    
    for k = 1:ndata
      freqindx{k,1} = [];
    end
  else
    error('incorrect specification of cfg.frequency');
  end
end % if cfg.frequency

if isfield(cfg, 'foilim')
  if ~ischar(cfg.foilim) && numel(cfg.foilim)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    minfreq = min(freqaxis);
    maxfreq = max(freqaxis);
    if all(cfg.foilim<minfreq) || all(cfg.foilim>maxfreq)
      error('the selected range falls outside the frequency axis in the data');
    end
    fbin = nan(1,2);
    fbin(1) = nearest(freqaxis, cfg.foilim(1), false, false);
    fbin(2) = nearest(freqaxis, cfg.foilim(2), false, false);
    cfg.foilim = freqaxis(fbin);
    
    for k = 1:ndata
      freqindx{k,1} = indx(fbin(1):fbin(2), k);
    end
    
  else
    error('incorrect specification of cfg.foilim');
  end
end % cfg.foilim

for k = 1:ndata
  if isequal(freqindx{k}, 1:length(varargin{k}.freq))
    % the cfg was updated, but no selection is needed for the data
    freqindx{k} = nan;
  end
end

end % function getselection_freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rptindx, cfg, rptdim, rptindxtap] = getselection_rpt(cfg, data, varargin)
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
      rptindxtap = find(tapers);
      [srt,ix]   = sort(tapers(tapers~=0));
      rptindxtap = rptindxtap(ix);
      %       cfg.trials = rptindx;
      % TODO FIXME think about whether this is a good or a bad thing...
      %warning('cfg.trials accounts for the number of tapers now');
    else
      rptindxtap = rptindx;
    end
    
    if ~isempty(rptindx) && rptindx(1)<1
      error('cannot select rpt/subj/rpttap smaller than 1');
    elseif ~isempty(rptindx) && rptindx(end)>rptsiz
      error('cannot select rpt/subj/rpttap larger than the number of repetitions in the data');
    end
    
    % commented out because of rpttap dilemma...
    %     cfg.trials = rptindx;
    
    return
  end
  
else
  rptindx = nan;
  rptindxtap = nan;
  
  % recover the rptdim if possible
  dimtok = tokenize(data.dimord, '_');
  rptdim = find(strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));
end % if isfield cfg.trials

end % function getselection_rpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posindx, cfg] = getselection_pos(cfg, data)
% possible specifications are <none>
posindx = 1:size(data.pos,1);

end % function getselection_pos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
