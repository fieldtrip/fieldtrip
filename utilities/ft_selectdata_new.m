function [varargout] = ft_selectdata_new(cfg, varargin)

% FT_SELECTDATA_NEW is deprecated, please use FT_SELECTDATA instead.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old documentation for reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% Use as
%  [data] = ft_selectdata_new(cfg, data, ...)
%
% The cfg artument is a configuration structure which can contain
%   cfg.tolerance   = scalar, tolerance value to determine equality of time/frequency bins (default = 1e-5)
%
% For data with trials or subjects as repetitions, you can specify
%   cfg.trials      = 1xN, trial indices to keep, can be 'all'. You can use logical indexing, where false(1,N) removes all the trials
%   cfg.avgoverrpt  = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a channel dimension you can specify
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%   cfg.avgoverchan = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a time dimension you can specify
%   cfg.latency     = scalar      -> can be 'all'
%   cfg.latency     = [beg end]
%   cfg.avgovertime = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a frequency dimension you can specify
%   cfg.frequency   = scalar    -> can be 'all'
%   cfg.frequency   = [beg end] -> this is less common, preferred is to use foilim
%   cfg.foilim      = [beg end]
%   cfg.avgoverfreq = string, can be 'yes' or 'no' (default = 'no')
%
% If multiple input arguments are provided, FT_SELECTDATA will adjust the
% individual inputs such that either the intersection across inputs is
% retained (i.e. only the channel/time/frequency points that are shared
% across all input arguments), or the union across inputs is retained
% (replacing missing data with nans). In either case, the order (e.g. of
% the channel labels) is made consistent across inputs. Multiple inputs in
% combination with the selection of trials is not supported. The exact
% behavior can be specified with
%   cfg.select      = 'intersect' or 'union' (default = 'intersect')

% Copyright (C) 2012-2014, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble provenance        % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar varargin  % this reads the input data in case the user specified the cfg.inputfile option

% determine the characteristics of the input data
dtype = ft_datatype(varargin{1});
for i=2:length(varargin)
  % ensure that all subsequent inputs are of the same type
  ok = ft_datatype(varargin{i}, dtype);
  if ~ok, ft_error('input data should be of the same datatype'); end
end

cfg = ft_checkconfig(cfg, 'renamed', {'selmode',  'select'});
cfg = ft_checkconfig(cfg, 'renamed', {'toilim' 'latency'});
cfg = ft_checkconfig(cfg, 'renamed', {'avgoverroi' 'avgoverpos'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'avg.pow' 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'avg.mom' 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'avg.nai' 'nai'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'trial.pow' 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'trial.mom' 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter' 'trial.nai' 'nai'});

cfg.tolerance = ft_getopt(cfg, 'tolerance', 1e-5);        % default tolerance for checking equality of time/freq axes
cfg.select    = ft_getopt(cfg, 'select',    'intersect'); % default is to take intersection, alternative 'union'
cfg.parameter = ft_getopt(cfg, 'parameter', {});

% this function only works for the upcoming (not yet standard) source representation without sub-structures
% update the old-style beamformer source reconstruction to the upcoming representation
if strcmp(dtype, 'source')
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
  end
end

if length(varargin)>1 && isfield(cfg, 'trials') && ~isequal(cfg.trials, 'all')
  ft_error('it is ambiguous to make a subselection of trials while at the same time concatenating multiple data structures')
end

if strcmp(cfg.select, 'union') && any(strcmp(dtype, {'raw', 'comp', 'source'}))
  ft_error('cfg.select ''union'' is not yet supported for %s data', dtype);
end

if ft_datatype(varargin{1}, 'raw')
  
  cfg.channel = ft_getopt(cfg, 'channel', 'all', 1); % empty definition by user is meaningful
  cfg.latency = ft_getopt(cfg, 'latency', 'all', 1);
  cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1);
  
  for i=1:length(varargin)
    varargin{i} = selfromraw(varargin{i}, 'rpt', cfg.trials, 'chan', cfg.channel, 'latency', cfg.latency);
  end
  
else % not raw or comp
  
  cfg.channel = ft_getopt(cfg, 'channel', 'all', 1);
  cfg.latency = ft_getopt(cfg, 'latency', 'all', 1);
  cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1);
  
  if ~isfield(cfg, 'foilim')
    cfg.frequency = ft_getopt(cfg, 'frequency', 'all', 1);
  end
  
  if isempty(cfg.parameter) && isfield(varargin{1}, 'dimord')
    dimord = varargin{1}.dimord;
  elseif ischar(cfg.parameter) && isfield(varargin{1}, [cfg.parameter 'dimord'])
    dimord = varargin{1}.([cfg.parameter 'dimord']);
  elseif ischar(cfg.parameter) && isfield(varargin{1}, 'dimord')
    dimord = varargin{1}.dimord;
  else
    ft_error('cannot determine which parameter to select from the data, please specify cfg.parameter');
  end
  
  dimtok = tokenize(dimord, '_');
  
  if isempty(cfg.parameter) || isequal(cfg.parameter ,'all')
    
    dimsiz = nan(size(dimtok));
    dimfields = cell(size(dimtok));
    
    % determine the size of each of the dimensions
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
        case '{pos}'
          dimsiz(i) = size(varargin{1}.pos,1);
          dimfields{i} = '{pos}';
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
        case 'ori'
          % the number of elements along this dimension is implicit
          dimsiz(i) = nan;
          dimfields{i} = 'implicit';
          
        case 'comp'
          ft_error('FIXME');
          
        case 'refchan'
          ft_error('FIXME');
          
        case 'voxel'
          ft_error('FIXME');
          
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
    
    % deal with the data dimensions whose size is only implicitly represented
    if any(strcmp(dimfields, 'implicit'))
      fn  = fieldnames(varargin{1})';
      for i=1:numel(fn)
        val = varargin{1}.(fn{i});
        siz = cellmatsize(val);
        clear val
        if isequalwithoutnans(siz, dimsiz)
          fprintf('using the "%s" field to determine the size along the unknown dimensions\n', fn{i});
          % update the size of all dimensions
          dimsiz = size(varargin{1}.(fn{i}));
          % update the fieldname of each dimension
          dimfields(strcmp(dimfields, 'implicit')) = dimtok(strcmp(dimfields, 'implicit'));
          break
        end
      end
      if any(strcmp(dimfields, 'implicit'))
        % it failed
        ft_error('could not determine the size of the implicit "%s" dimension', dimfields{strcmp(dimfields, 'implicit')});
      end
    end
    
    % select the fields based on the dimord
    fn  = fieldnames(varargin{1})'; % it should be a row-array
    fn  = setdiff(fn, {'pos', 'label', 'time', 'freq', 'cfg', 'hdr', 'grad', 'elec'});
    sel = false(size(fn));
    for i=1:numel(fn)
      sel(i) = isequal(size(varargin{1}.(fn{i})), dimsiz) || isequal(size(varargin{1}.(fn{i})), [dimsiz 1]);
    end
    cfg.parameter = fn(sel);
    
    clear dimsiz dimfields
    
  end % is isempty(cfg.parameter)
  
  % these are the fields in which the selection will be made
  datfields = cfg.parameter;
  if ~iscell(datfields)
    datfields = {datfields};
  end
  
  hasrpt     = any(ismember(dimtok, {'rpt', 'subj'}));
  hasrpttap  = any(ismember(dimtok, 'rpttap'));
  haspos     = any(ismember(dimtok, {'pos', '{pos}'}));
  haschan    = any(ismember(dimtok, 'chan'));
  haschancmb = any(ismember(dimtok, 'chancmb'));
  hasfreq    = any(ismember(dimtok, 'freq'));
  hastime    = any(ismember(dimtok, 'time'));
  
  haspos     = haspos     && isfield(varargin{1}, 'pos');
  haschan    = haschan    && isfield(varargin{1}, 'label');
  haschancmb = haschancmb && isfield(varargin{1}, 'labelcmb');
  hasfreq    = hasfreq    && isfield(varargin{1}, 'freq');
  hastime    = hastime    && isfield(varargin{1}, 'time');
  
  avgoverpos  = istrue(ft_getopt(cfg, 'avgoverpos',  false)); % at some places it is also referred to as roi (region-of-interest)
  avgoverrpt  = istrue(ft_getopt(cfg, 'avgoverrpt',  false));
  avgoverchan = istrue(ft_getopt(cfg, 'avgoverchan', false));
  avgoverfreq = istrue(ft_getopt(cfg, 'avgoverfreq', false));
  avgovertime = istrue(ft_getopt(cfg, 'avgovertime', false));
  
  if avgoverpos,  assert(haspos,  'there are no source positions, so averaging is not possible'); end
  if avgoverrpt,  assert(hasrpt||hasrpttap, 'there are no repetitions, so averaging is not possible'); end
  if avgoverchan, assert(haschan, 'there is no channel dimension, so averaging is not possible'); end
  if avgoverfreq, assert(hasfreq, 'there is no frequency dimension, so averaging is not possible'); end
  if avgovertime, assert(hastime, 'there is no time dimension, so averaging over time is not possible'); end
  
  % by default we keep most of the dimensions in the data structure when averaging over them
  keeprptdim  = istrue(ft_getopt(cfg, 'keeprptdim', false));
  keepposdim  = istrue(ft_getopt(cfg, 'keepposdim',  true));
  keepchandim = istrue(ft_getopt(cfg, 'keepchandim', true));
  keepfreqdim = istrue(ft_getopt(cfg, 'keepfreqdim', true));
  keeptimedim = istrue(ft_getopt(cfg, 'keeptimedim', true));
  
  if strcmp(cfg.select, 'union') && (avgoverpos || avgoverrpt || avgoverchan || avgoverfreq || avgovertime)
    ft_error('cfg.select ''union'' in combination with averaging across one of the dimensions is not implemented');
  end
  
  if avgoverpos
    for i=1:length(varargin)
      % must be a source representation, not a volume representation
      varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source');
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 2:
  %   ensure that the cfg is fully contained in the data and consistent over all inputs
  %   get the selection along each of the dimensions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % FIXMEroboos this implementation is not yet complete
  % dtype = 'new';
  
  switch dtype
    % this switch-list is consistent with ft_datatype
    
    case {'new'}
      % FIXMEroboos this implementation is not yet complete
      % trim the selection to all inputs
      if haspos,     [selpos,     cfg]  = getselection_pos   (cfg, varargin{:}, cfg.tolerance, cfg.select); end
      if haschan,    [selchan,    cfg] = getselection_chan   (cfg, varargin{:}, cfg.select); end
      if haschancmb, [selchancmb, cfg] = getselection_chancmb(cfg, varargin{:}, cfg.select); end
      if hastime,    [seltime,    cfg] = getselection_time   (cfg, varargin{:}, cfg.tolerance, cfg.select); end
      if hasfreq,    [selfreq,    cfg] = getselection_freq   (cfg, varargin{:}, cfg.tolerance, cfg.select); end
      
      
      for i=1:numel(varargin)
        % the rpt selection should only work with a single data argument
        % in case tapers were kept, selrpt~=selrpttap, otherwise selrpt==selrpttap
        [selrpt{i}, dum, rptdim{i}, selrpttap{i}] = getselection_rpt(cfg, varargin{i}, 'datfields', datfields);
        
        if haspos,     varargin{i} = makeselection(varargin{i}, find(ismember(dimtok, {'pos', '{pos}'})), selpos{i},     avgoverpos,  datfields, cfg.select); end
        if haschan,    varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')),              selchan{i},    avgoverchan, datfields, cfg.select); end
        if haschancmb, varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chancmb')),           selchancmb{i}, false,       datfields, cfg.select); end
        if hastime,    varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')),              seltime{i},    avgovertime, datfields, cfg.select); end
        if hasfreq,    varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'freq')),              selfreq{i},    avgoverfreq, datfields, cfg.select); end
        if hasrpt,     varargin{i} = makeselection(varargin{i}, find(ismember(dimtok,{'rpt', 'subj'})),   selrpt{i},     avgoverrpt,  datfields, 'intersect'); end
        if hasrpttap,  varargin{i} = makeselection(varargin{i}, rptdim{i},                                selrpttap{i},  avgoverrpt,  datfields, 'intersect'); end
        
        if haspos,  varargin{i} = makeselection_pos(varargin{i},  selpos{i},  avgoverpos);  end % update the pos field
        if haschan, varargin{i} = makeselection_chan(varargin{i}, selchan{i}, avgoverchan); end % update the label field
        if hastime, varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); end % update the time field
        if hasfreq, varargin{i} = makeselection_freq(varargin{i}, selfreq{i}, avgoverfreq); end % update the time field
        if hasrpt || hasrpttap, varargin{i} = makeselection_rpt (varargin{i}, selrpt{i});   end % avgoverrpt for the supporting fields is dealt with later
        
        % also deal with the supporting cumtapcnt field, because it has a frequency dimension when time dimension is present
        % this is a temporary workaround, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2509
        if isfield(varargin{i}, 'cumtapcnt') && hastime
          varargin{i} = makeselection_cumtapcnt(varargin{i}, selfreq{i}, avgoverfreq);
        end
        
        % make an exception for the covariance here (JM 20131128)
        if isfield(varargin{i}, 'cov') && (all(~isnan(selrpt{i})) || all(~isnan(selchan{i})))
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'chan'))+[0 1], selchan{i}, avgoverchan, {'cov'}, cfg.select);
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'rpt')),        selrpt{i},  avgoverrpt,  {'cov'}, 'intersect');
          datfields   = [datfields {'cov'}];
        end
        
      end % for varargin
      
      % in the case of selmode='union', create the union of the descriptive axes
      if strcmp(cfg.select, 'union')
        
        if haschan
          label = varargin{1}.label;
          for i=2:numel(varargin)
            tmplabel   = varargin{i}.label;
            emptylabel = find(cellfun('isempty', label));
            for k=emptylabel(:)'
              label{k} = tmplabel{k};
            end
          end
          for i=1:numel(varargin)
            varargin{i}.label = label;
          end
        end % haschan
        
        if hastime
          time = varargin{1}.time;
          for i=2:numel(varargin)
            tmptime = varargin{i}.time;
            time(~isfinite(time)) = tmptime(~isfinite(time));
          end
          for i=1:numel(varargin)
            varargin{i}.time  = time;
          end
        end % hastime
        
      end % select=union
      
    case 'timelock'
      % trim the selection to all inputs
      [selchan, cfg] = getselection_chan(cfg, varargin{:}, cfg.select);
      [seltime, cfg] = getselection_time(cfg, varargin{:}, cfg.tolerance, cfg.select);
      
      selrpt = cell(numel(varargin),1);
      for i=1:numel(varargin)
        [selrpt{i}] = getselection_rpt (cfg, varargin{i}, 'datfields', datfields);
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')),                      selchan{i}, avgoverchan, datfields, cfg.select);
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')),                      seltime{i}, avgovertime, datfields, cfg.select);
        varargin{i} = makeselection(varargin{i}, find(ismember(dimtok,{'rpt', 'rpttap', 'subj'})), selrpt{i},  avgoverrpt,  datfields, 'intersect');
        
        varargin{i} = makeselection_chan(varargin{i}, selchan{i}, avgoverchan); % update the label field
        varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); % update the time field
        varargin{i} = makeselection_rpt (varargin{i}, selrpt{i}); % avgoverrpt for the supporting fields is dealt with later
        
        % make an exception for the covariance here (JM 20131128)
        if isfield(varargin{i}, 'cov') && (all(~isnan(selrpt{i})) || all(~isnan(selchan{i})))
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'chan'))+[0 1], selchan{i}, avgoverchan, {'cov'}, cfg.select);
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok, 'rpt')),        selrpt{i},  avgoverrpt,  {'cov'}, 'intersect');
          datfields   = [datfields {'cov'}];
        end
        
      end % varargin
      
      % in the case of selmode='union', create the union of the descriptive axes
      if strcmp(cfg.select, 'union')
        label = varargin{1}.label;
        time  = varargin{1}.time;
        
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
      [selchan, cfg] = getselection_chan(cfg, varargin{:},                cfg.select); % tolerance not needed
      [selfreq, cfg] = getselection_freq(cfg, varargin{:}, cfg.tolerance, cfg.select); % freq is always present
      if hastime, [seltime, cfg] = getselection_time(cfg, varargin{:}, cfg.tolerance, cfg.select); end
      
      selrpt    = cell(numel(varargin),1);
      selrpttap = cell(numel(varargin),1);
      rptdim    = cell(numel(varargin),1);
      for i=1:numel(varargin)
        % the rpt selection stays within this loop, it only should work with a single data argument anyway
        % in case tapers were kept, selrpt~=selrpttap, otherwise selrpt==selrpttap
        [selrpt{i}, dum, rptdim{i}, selrpttap{i}] = getselection_rpt(cfg, varargin{i}, 'datfields', datfields);
        
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'chan')), selchan{i},   avgoverchan, datfields, cfg.select);
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'freq')), selfreq{i},   avgoverfreq, datfields, cfg.select);
        varargin{i} = makeselection(varargin{i}, rptdim{i},                   selrpttap{i}, avgoverrpt,  datfields, 'intersect');
        
        varargin{i} = makeselection_chan(varargin{i}, selchan{i}, avgoverchan); % update the label field
        varargin{i} = makeselection_freq(varargin{i}, selfreq{i}, avgoverfreq); % update the freq field
        varargin{i} = makeselection_rpt(varargin{i},  selrpt{i}); % avgoverrpt for the supporting fields is dealt with later
        
        if hastime
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime{i}, avgovertime, datfields, cfg.select);
          varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); % update the time field
        end
        
        % also deal with the supporting cumtapcnt field, because it has a frequency dimension when time dimension is present
        % this is a temporary workaround, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2509
        if hastime && isfield(varargin{i}, 'cumtapcnt')
          varargin{i} = makeselection_cumtapcnt(varargin{i}, selfreq{i}, avgoverfreq);
        end
      end % varargin
      
      % in the case of selmode='union', create the union of the descriptive axes
      if strcmp(cfg.select, 'union')
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
      % trim the selection to all inputs
      [selpos, cfg] = getselection_pos(cfg, varargin{:}, cfg.tolerance, cfg.select);
      if hastime, [seltime, cfg] = getselection_time(cfg, varargin{:}, cfg.tolerance, cfg.select); end
      if hasfreq, [selfreq, cfg] = getselection_freq(cfg, varargin{:}, cfg.tolerance, cfg.select); end
      
      for i=1:numel(varargin)
        % get the selection from all inputs
        varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'pos') | strcmp(dimtok,'{pos}')), selpos{i}, avgoverpos, datfields, cfg.select);
        varargin{i} = makeselection_pos(varargin{i}, selpos{i}, avgoverpos); % update the pos field
        
        % FIXME this code does not deal with repetitions
        
        if hastime,
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'time')), seltime{i}, avgovertime, datfields, cfg.select);
          varargin{i} = makeselection_time(varargin{i}, seltime{i}, avgovertime); % update the time field
        end
        if hasfreq,
          varargin{i} = makeselection(varargin{i}, find(strcmp(dimtok,'freq')), selfreq{i},   avgoverfreq, datfields, cfg.select);
          varargin{i} = makeselection_freq(varargin{i}, selfreq{i}, avgoverfreq); % update the freq field
        end
        
      end % varargin
      
    case 'freqmvar'
      ft_error('FIXME');
      
    case 'mvar'
      ft_error('FIXME');
      
    case 'spike'
      ft_error('FIXME');
      
    case 'volume'
      ft_error('FIXME');
      
    case 'dip'
      ft_error('FIXME');
      
    case 'chan'
      % this results from avgovertime/avgoverfreq after timelockstatistics or freqstatistics
      ft_error('FIXME');
      
    otherwise
      % try to get the selection based on the field name
      seldim = cell(size(dimtok));
      for j=1:numel(seldim)
        seldim(j) = feval(['getselection_' dimtok{j}], cfg, varargin{i});
      end
  end % switch dtype
  
  % update the fields and the dimord
  keepdim   = true(size(dimtok));
  keepfield = unique(dimtok);
  sel = strcmp(keepfield, '{pos}'); if any(sel), keepfield(sel) = {'pos'}; end
  sel = strcmp(keepfield, 'chan');  if any(sel), keepfield(sel) = {'label'}; end
  
  if avgoverchan && ~keepchandim
    keepdim(strcmp(dimtok, 'chan')) = false;
    keepfield = setdiff(keepfield, 'label');
  else
    keepfield = [keepfield 'label'];
  end
  
  if avgoverfreq && ~keepfreqdim
    keepdim(strcmp(dimtok, 'freq')) = false;
    keepfield = setdiff(keepfield, 'freq');
  else
    keepfield = [keepfield 'freq'];
  end
  
  if avgovertime && ~keeptimedim
    keepdim(strcmp(dimtok, 'time')) = false;
    keepfield = setdiff(keepfield, 'time');
  else
    keepfield = [keepfield 'time'];
  end
  
  if avgoverpos && ~keepposdim
    keepdim(strcmp(dimtok, 'pos'))   = false;
    keepdim(strcmp(dimtok, '{pos}')) = false;
    keepfield = setdiff(keepfield, {'pos' '{pos}' 'dim'});
  else
    keepfield = [keepfield {'pos' '{pos}' 'dim'}];
  end
  
  if avgoverrpt && ~keeprptdim
    keepdim(ismember(dimtok, {'rpt', 'rpttap', 'subj'})) = false;
    keepfield = setdiff(keepfield, {'cumtapcnt' 'cumsumcnt' 'sampleinfo' 'trialinfo'});
  else
    keepfield = [keepfield {'cumtapcnt' 'cumsumcnt' 'sampleinfo' 'trialinfo'}];
  end
  
  % remove all fields from the dimord that do not pertain to the selection
  for i=1:numel(varargin)
    varargin{i}.dimord = sprintf('%s_', dimtok{keepdim});
    varargin{i}.dimord = varargin{i}.dimord(1:end-1);  % remove the last '_'
  end
  
  for i=1:numel(varargin)
    for j=1:numel(datfields)
      varargin{i}.(datfields{j}) = squeezedim(varargin{i}.(datfields{j}), ~keepdim);
    end
  end
  
  % remove all fields from the data that do not pertain to the selection
  for i=1:numel(varargin)
    varargin{i} = keepfields(varargin{i}, [datfields {'cfg' 'dimord' 'elec' 'grad'} keepfield]);
  end
  
end % if raw or something else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3:
%   if desired, concatenate over repetitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout = varargin;

ft_postamble debug              % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig        % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance         % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
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

end % main function ft_selectdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
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
      % the selindx value of NaN indicates that it is not needed to make a selection
      if isempty(selindx) || all(~isnan(selindx))
        data.(datfields{i}) = cellmatselect(data.(datfields{i}), seldim, selindx);
      end
      if avgoverdim
        data.(datfields{i}) = cellmatmean(data.(datfields{i}), seldim);
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
          ft_error('unsupported dimension (%d) for making a selection for %s', seldim, datfields{i});
      end
    end
    if avgoverdim
      data.(datfields{i}) = mean(data.(datfields{i}), seldim);
    end
end % switch

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
  % compute the mean frequency
  if ~isnan(selfreq)
    data.freq = mean(data.freq(selfreq));
  else
    data.freq = mean(data.freq);
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
function data = makeselection_time(data, seltime, avgovertime)
if avgovertime
  % compute the mean latency
  if ~isnan(seltime)
    data.time = mean(data.time(seltime));
  else
    data.time = mean(data.time);
  end
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
function data = makeselection_cumtapcnt(data, selfreq, avgoverfreq)

if ~isfield(data, 'time')
  ft_error('the subfunction makeselection_cumtapcnt should only be called when there is a time dimension in the data');
end
if ~isfield(data, 'cumtapcnt')
  return;
end

if avgoverfreq
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
function data = makeselection_pos(data, selpos, avgoverpos)
if ~isnan(selpos)
  data.pos = data.pos(selpos, :);
end
if avgoverpos
  data.pos = mean(data.pos, 1);
end
end % function makeselection_pos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chanindx, cfg] = getselection_chan(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

% loop over data once to initialize
chanindx = cell(ndata,1);
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
      ft_error('invalid value for cfg.select');
  end % switch
  
  for k = 1:ndata
    chanindx{k,1} = indx(:,k);
  end
  cfg.channel = label;
  
else
  for k = 1:ndata
    % the nan return value specifies that no selection was specified
    chanindx{k,1} = nan;
  end
  
end

end % function getselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chancmbindx, cfg] = getselection_chancmb(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

chancmbindx = cell(ndata,1);

if isfield(cfg, 'channelcmb')
  for k = 1:ndata
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb, varargin{k}.labelcmb);
  end
  
  ft_error('selection of channelcmb is not yet implemented');
  
else
  for k = 1:ndata
    % the nan return value specifies that no selection was specified
    chancmbindx{k} = nan;
  end
end
end % function getselection_chancmb

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
  
  % the nan return value specifies that no selection was specified
  timeindx{k,1} = nan;
  
  % update the axis along which the frequencies are defined
  timeaxis = union(timeaxis, round(varargin{k}.time(:)/tol)*tol);
end

indx = nan+zeros(numel(timeaxis), ndata);
for k = 1:ndata
  [dum, ix, iy] = intersect(timeaxis, round(varargin{k}.time(:)/tol)*tol);
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
    ft_error('invalid value for cfg.select');
end

if isfield(cfg, 'latency')
  % deal with string selection
  if ischar(cfg.latency)
    if strcmp(cfg.latency, 'all')
      cfg.latency = [min(timeaxis) max(timeaxis)];
    else
      ft_error('incorrect specification of cfg.latency');
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
      ft_error('the selected time range falls outside the time axis in the data');
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
    ft_error('incorrect specification of cfg.latency');
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
  
  % the nan return value specifies that no selection was specified
  freqindx{k,1} = nan;
  
  % update the axis along which the frequencies are defined
  freqaxis = union(freqaxis, round(varargin{k}.freq(:)/tol)*tol);
end

indx = nan+zeros(numel(freqaxis), ndata);
for k = 1:ndata
  [dum, ix, iy] = intersect(freqaxis, round(varargin{k}.freq(:)/tol)*tol);
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
    ft_error('invalid value for cfg.select');
end

if isfield(cfg, 'frequency')
  % deal with string selection
  if ischar(cfg.frequency)
    if strcmp(cfg.frequency, 'all')
      cfg.frequency = [min(freqaxis) max(freqaxis)];
    else
      ft_error('incorrect specification of cfg.frequency');
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
      ft_error('the selected range falls outside the frequency axis in the data');
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
    ft_error('incorrect specification of cfg.frequency');
  end
end % if cfg.frequency

if isfield(cfg, 'foilim')
  if ~ischar(cfg.foilim) && numel(cfg.foilim)==2
    % the [min max] range can be specifed with +inf or -inf, but should
    % at least partially overlap with the time axis of the input data
    minfreq = min(freqaxis);
    maxfreq = max(freqaxis);
    if all(cfg.foilim<minfreq) || all(cfg.foilim>maxfreq)
      ft_error('the selected range falls outside the frequency axis in the data');
    end
    fbin = nan(1,2);
    fbin(1) = nearest(freqaxis, cfg.foilim(1), false, false);
    fbin(2) = nearest(freqaxis, cfg.foilim(2), false, false);
    cfg.foilim = freqaxis(fbin);
    
    for k = 1:ndata
      freqindx{k,1} = indx(fbin(1):fbin(2), k);
    end
    
  else
    ft_error('incorrect specification of cfg.foilim');
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

% start with the initual guess for the dimord
if isfield(data, 'dimord')
  dimord = data.dimord;
end

% perhaps there is a specific dimord for the data fields of interest
for i=1:length(datfields)
  if isfield(data, [datfields{i} 'dimord'])
    dimord = data.([datfields{i} 'dimord']);
    break
  end
end

dimtok = tokenize(dimord, '_');

if isfield(cfg, 'trials') && ~isequal(cfg.trials, 'all') && ~isempty(datfields)
  
  rptdim = find(strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));
  rptindx    = nan; % the nan return value specifies that no selection was specified
  rptindxtap = nan; % the nan return value specifies that no selection was specified
  
  if isempty(rptdim)
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
      % account for the tapers
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
      ft_error('cannot select rpt/subj/rpttap smaller than 1');
    elseif ~isempty(rptindx) && rptindx(end)>rptsiz
      ft_error('cannot select rpt/subj/rpttap larger than the number of repetitions in the data');
    end
    
    % commented out because of rpttap dilemma...
    %     cfg.trials = rptindx;
    
    return
  end
  
else
  % recover the rptdim if possible
  rptdim     = find(strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));
  rptindx    = nan;
  rptindxtap = nan;
  
end % if isfield cfg.trials

end % function getselection_rpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posindx, cfg] = getselection_pos(cfg, varargin)
% possible specifications are <none>

ndata = numel(varargin)-2;
tol   = varargin{end-1}; % FIXME this is still ignored
selmode = varargin{end}; % FIXME this is still ignored

for i=2:ndata
  if ~isequal(varargin{i}.pos, varargin{1}.pos)
    % FIXME it would be possible here to make a selection based on intersect or union
    ft_error('source positions are different');
  end
end % for
for i=1:ndata
  posindx{i} = nan;    % the nan return value specifies that no selection was specified
end
end % function getselection_pos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = squeezedim(x, dim)
siz = size(x);
for i=(numel(siz)+1):numel(dim)
  % all trailing singleton dimensions have length 1
  siz(i) = 1;
end
x = reshape(x, [siz(~dim) 1]);
end % function squeezedim

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the size of data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function siz = cellmatsize(x)
if iscell(x)
  cellsize = numel(x);          % the number of elements in the cell-array
  [dum, indx] = max(cellfun(@numel, x));
  matsize = size(x{indx});      % the size of the content of the cell-array
  siz  = [cellsize matsize];    % concatenate the two
else
  siz = size(x);
end
end % function cellmatsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to make a selextion in data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cellmatselect(x, seldim, selindx)
if iscell(x)
  if seldim==1
    x = x(selindx);
  else
    for i=1:numel(x)
      switch seldim
        case 2
          x{i} = x{i}(selindx,:,:,:,:);
        case 3
          x{i} = x{i}(:,selindx,:,:,:);
        case 4
          x{i} = x{i}(:,:,selindx,:,:);
        case 5
          x{i} = x{i}(:,:,:,selindx,:);
        case 6
          x{i} = x{i}(:,:,:,:,selindx);
        otherwise
          ft_error('unsupported dimension (%d) for making a selection', seldim);
      end % switch
    end % for
  end
else
  switch seldim
    case 1
      x = x(selindx,:,:,:,:,:);
    case 2
      x = x(:,selindx,:,:,:,:);
    case 3
      x = x(:,:,selindx,:,:,:);
    case 4
      x = x(:,:,:,selindx,:,:);
    case 5
      x = x(:,:,:,:,selindx,:);
    case 6
      x = x(:,:,:,:,:,selindx);
    otherwise
      ft_error('unsupported dimension (%d) for making a selection', seldim);
  end
end
end % function cellmatselect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to take an average in data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cellmatmean(x, seldim)
if iscell(x)
  if seldim==1
    for i=2:numel(x)
      x{1} = x{1} + x{i};
    end
    x = {x{1}/numel(x)};
  else
    for i=1:numel(x)
      x{i} = mean(x{i}, seldim-1);
    end % for
  end
else
  x = mean(x, seldim);
end
end % function cellmatmean

function dimord = paramdimord(data, param)
if isfield(data, [param 'dimord'])
  dimord = data.([param 'dimord']);
else
  dimord = data.dimord;
end
end % function paramdimord

function dimtok = paramdimtok(data, param)
dimord = paramdimord(data, param);
dimtok = tokenize(dimord, '_');
end % function paramdimtok

