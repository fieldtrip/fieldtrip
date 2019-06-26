function [varargout] = ft_selectdata(cfg, varargin)

% FT_SELECTDATA makes a selection in the input data along specific data
% dimensions, such as channels, time, frequency, trials, etc. It can also
% be used to average the data along each of the specific dimensions.
%
% Use as
%  [data] = ft_selectdata(cfg, data, ...)
%
% The cfg argument is a configuration structure which can contain
%   cfg.tolerance   = scalar, tolerance value to determine equality of time/frequency bins (default = 1e-5)
%
% For data with trials or subjects as repetitions, you can specify
%   cfg.trials      = 1xN, trial indices to keep, can be 'all'. You can use logical indexing, where false(1,N) removes all the trials
%   cfg.avgoverrpt  = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a channel dimension you can specify
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%   cfg.avgoverchan = string, can be 'yes' or 'no' (default = 'no')
%   cfg.nanmean     = string, can be 'yes' or 'no' (default = 'no')
%
% For data with channel combinations you can specify
%   cfg.channelcmb     = Nx2 cell-array with selection of channels (default = 'all'), see FT_CHANNELCOMBINATION
%   cfg.avgoverchancmb = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a time dimension you can specify
%   cfg.latency     = scalar or string, can be 'all', 'minperiod', 'maxperiod', 'prestim', 'poststim', or [beg end], specify time range in seconds
%   cfg.avgovertime = string, can be 'yes' or 'no' (default = 'no')
%   cfg.nanmean     = string, can be 'yes' or 'no' (default = 'no')
%
% For data with a frequency dimension you can specify
%   cfg.frequency   = scalar or string, can be 'all', or [beg end], specify frequency range in Hz
%   cfg.avgoverfreq = string, can be 'yes' or 'no' (default = 'no')
%   cfg.nanmean     = string, can be 'yes' or 'no' (default = 'no')
%
% If multiple input arguments are provided, FT_SELECTDATA will adjust the individual inputs
% such that either the intersection across inputs is retained (i.e. only the channel, time,
% and frequency points that are shared across all input arguments), or that the union across
% inputs is retained (replacing missing data with nans). In either case, the order (e.g. of
% the channels) is made consistent across inputs.  The behavior can be specified with
%   cfg.select      = string, can be 'intersect' or 'union' (default = 'intersect')
%
% See also FT_DATATYPE, FT_CHANNELSELECTION, FT_CHANNELCOMBINATION

% Undocumented options
%   cfg.keeprptdim     = 'yes' or 'no'
%   cfg.keepposdim     = 'yes' or 'no'
%   cfg.keepchandim    = 'yes' or 'no'
%   cfg.keepchancmbdim = 'yes' or 'no'
%   cfg.keepfreqdim    = 'yes' or 'no'
%   cfg.keeptimedim    = 'yes' or 'no'
%   cfg.avgoverpos

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

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble trackconfig
ft_preamble loadvar varargin
ft_preamble provenance varargin

% determine the characteristics of the input data
dtype = ft_datatype(varargin{1});
for i=2:length(varargin)
  % ensure that all subsequent inputs are of the same type
  ok = ft_datatype(varargin{i}, dtype);
  if ~ok, ft_error('input data should be of the same datatype'); end
end

% this only works with certain data types, it is not meant for descriptive fields such as elec, grad, opto, layout, etc.
assert(~ismember(dtype, {'elec', 'grad', 'opto', 'layout'}), 'invalid input data type "%s"', dtype);

% ensure that the user does not give invalid selection options
cfg = ft_checkconfig(cfg, 'forbidden', {'foi', 'toi'});

cfg = ft_checkconfig(cfg, 'renamed',    {'selmode',    'select'});
cfg = ft_checkconfig(cfg, 'renamed',    {'toilim',     'latency'});
cfg = ft_checkconfig(cfg, 'renamed',    {'foilim',     'frequency'});
cfg = ft_checkconfig(cfg, 'renamed',    {'avgoverroi', 'avgoverpos'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'avg.nai', 'nai'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'trial.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'trial.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'parameter', 'trial.nai', 'nai'});

cfg.tolerance = ft_getopt(cfg, 'tolerance', 1e-5);        % default tolerance for checking equality of time/freq axes
cfg.select    = ft_getopt(cfg, 'select',   'intersect');  % default is to take intersection, alternative 'union'

if strcmp(dtype, 'volume') || strcmp(dtype, 'segmentation')
  % it must be a source representation, not a volume representation
  for i=1:length(varargin)
    varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source');
  end
  dtype = 'source';
else
  % check that the data is according to the latest FieldTrip representation
  for i=1:length(varargin)
    varargin{i} = ft_checkdata(varargin{i});
  end
end

% this function only works for the upcoming (not yet standard) source representation without sub-structures
% update the old-style beamformer source reconstruction to the upcoming representation
if strcmp(dtype, 'source')
  if isfield(varargin{1}, 'avg') && isstruct(varargin{1}.avg)
    restoreavg = fieldnames(varargin{1}.avg);
  else
    restoreavg = {};
  end
  for i=1:length(varargin)
    varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
  end
end

cfg.latency = ft_getopt(cfg, 'latency', 'all', 1);
if isnumeric(cfg.latency) && numel(cfg.latency)==2 && cfg.latency(1)==cfg.latency(2)
  % this is better specified by a single number
  cfg.latency = cfg.latency(1);
end

cfg.channel = ft_getopt(cfg, 'channel', 'all', 1);
cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1);

if length(varargin)>1 && ~isequal(cfg.trials, 'all')
  ft_error('it is ambiguous to make a subselection of trials while at the same time concatenating multiple data structures')
end

cfg.frequency = ft_getopt(cfg, 'frequency', 'all', 1);
if isnumeric(cfg.frequency) && numel(cfg.frequency)==2 && cfg.frequency(1)==cfg.frequency(2)
  % this is better specified by a single number
  cfg.frequency = cfg.frequency(1);
end

datfield  = fieldnames(varargin{1});
for i=2:length(varargin)
  % only consider fields that are present in all inputs
  datfield = intersect(datfield, fieldnames(varargin{i}));
end
datfield  = setdiff(datfield, {'label' 'labelcmb'});  % these fields will be used for selection, but are not treated as data fields
datfield  = setdiff(datfield, {'dim'});               % not used for selection, also not treated as data field
datfield  = setdiff(datfield, ignorefields('selectdata'));
orgdim1   = datfield(~cellfun(@isempty, regexp(datfield, 'label$')) & cellfun(@isempty, regexp(datfield, '^csd'))); % xxxlabel, with the exception of csdlabel
datfield  = setdiff(datfield, orgdim1);
datfield  = datfield(:)';
orgdim1   = datfield(~cellfun(@isempty, regexp(datfield, 'dimord$'))); % xxxdimord
datfield  = setdiff(datfield, orgdim1);
datfield  = datfield(:)';

sel = strcmp(datfield, 'cumtapcnt');
if any(sel)
  % move this field to the end, as it is needed to make the selections in the other fields
  datfield(sel) = [];
  datfield = [datfield {'cumtapcnt'}];
end

orgdim2 = cell(size(orgdim1));
for i=1:length(orgdim1)
  orgdim2{i} = varargin{1}.(orgdim1{i});
end

dimord = cell(size(datfield));
for i=1:length(datfield)
  dimord{i} = getdimord(varargin{1}, datfield{i});
end

% do not consider fields of which the dimensions are unknown
% sel = cellfun(@isempty, regexp(dimord, 'unknown'));
% for i=find(~sel)
%   fprintf('not including "%s" in selection\n', datfield{i});
% end
% datfield = datfield(sel);
% dimord   = dimord(sel);

% determine all dimensions that are present in all data fields
dimtok = {};
for i=1:length(datfield)
  dimtok = cat(2, dimtok, tokenize(dimord{i}, '_'));
end
dimtok = unique(dimtok);

hasspike   = any(ismember(dimtok, 'spike'));
haspos     = any(ismember(dimtok, {'pos' '{pos}'}));
haschan    = any(ismember(dimtok, {'chan' '{chan}'}));
haschancmb = any(ismember(dimtok, 'chancmb'));
hasfreq    = any(ismember(dimtok, 'freq'));
hastime    = any(ismember(dimtok, 'time'));
hasrpt     = any(ismember(dimtok, {'rpt' 'subj' '{rpt}'}));
hasrpttap  = any(ismember(dimtok, 'rpttap'));

if hasspike
  % cfg.latency is used to select individual spikes, not to select from a continuously sampled time axis
  hastime = false;
end

clear dimtok

haspos     = haspos     && isfield(varargin{1}, 'pos');
haschan    = haschan    && isfield(varargin{1}, 'label');
haschancmb = haschancmb && isfield(varargin{1}, 'labelcmb');
hasfreq    = hasfreq    && isfield(varargin{1}, 'freq');
hastime    = hastime    && isfield(varargin{1}, 'time');

% do a sanity check on all input arguments
if haspos,     assert(all(cellfun(@isfield, varargin, repmat({'pos'},      size(varargin)))), 'not all input arguments have a "pos" field'); end
if haschan,    assert(all(cellfun(@isfield, varargin, repmat({'label'},    size(varargin)))), 'not all input arguments have a "label" field'); end
if haschancmb, assert(all(cellfun(@isfield, varargin, repmat({'labelcmb'}, size(varargin)))), 'not all input arguments have a "labelcmb" field'); end
if hasfreq,    assert(all(cellfun(@isfield, varargin, repmat({'freq'},     size(varargin)))), 'not all input arguments have a "freq" field'); end
if hastime,    assert(all(cellfun(@isfield, varargin, repmat({'time'},     size(varargin)))), 'not all input arguments have a "time" field'); end

avgoverpos     = istrue(ft_getopt(cfg, 'avgoverpos',  false)); % at some places it is also referred to as roi (region-of-interest)
avgoverchan    = istrue(ft_getopt(cfg, 'avgoverchan', false));
avgoverchancmb = istrue(ft_getopt(cfg, 'avgoverchancmb', false));
avgoverfreq    = istrue(ft_getopt(cfg, 'avgoverfreq', false));
avgovertime    = istrue(ft_getopt(cfg, 'avgovertime', false));
avgoverrpt     = istrue(ft_getopt(cfg, 'avgoverrpt',  false));

% do a sanity check for the averaging options
if avgoverpos,     assert(haspos,     'there are no source positions, so averaging is not possible'); end
if avgoverchan,    assert(haschan,    'there is no channel dimension, so averaging is not possible'); end
if avgoverchancmb, assert(haschancmb, 'there are no channel combinations, so averaging is not possible'); end
if avgoverfreq,    assert(hasfreq,    'there is no frequency dimension, so averaging is not possible'); end
if avgovertime,    assert(hastime,    'there is no time dimension, so averaging over time is not possible'); end
if avgoverrpt,     assert(hasrpt||hasrpttap, 'there are no repetitions, so averaging is not possible'); end

% set averaging function
cfg.nanmean = ft_getopt(cfg, 'nanmean', 'no');
if strcmp(cfg.nanmean, 'yes')
  average = @nanmean;
else
  average = @mean;
end

% by default we keep most of the dimensions in the data structure when averaging over them
keepposdim     = istrue(ft_getopt(cfg, 'keepposdim',  true));
keepchandim    = istrue(ft_getopt(cfg, 'keepchandim', true));
keepchancmbdim = istrue(ft_getopt(cfg, 'keepchancmbdim', true));
keepfreqdim    = istrue(ft_getopt(cfg, 'keepfreqdim', true));
keeptimedim    = istrue(ft_getopt(cfg, 'keeptimedim', true));
keeprptdim     = istrue(ft_getopt(cfg, 'keeprptdim', ~avgoverrpt));

if ~keepposdim,     assert(avgoverpos,     'removing a dimension is only possible when averaging'); end
if ~keepchandim,    assert(avgoverchan,    'removing a dimension is only possible when averaging'); end
if ~keepchancmbdim, assert(avgoverchancmb, 'removing a dimension is only possible when averaging'); end
if ~keepfreqdim,    assert(avgoverfreq,    'removing a dimension is only possible when averaging'); end
if ~keeptimedim,    assert(avgovertime,    'removing a dimension is only possible when averaging'); end
if ~keeprptdim,     assert(avgoverrpt,     'removing a dimension is only possible when averaging'); end

if strcmp(cfg.select, 'union') && (avgoverpos || avgoverchan || avgoverchancmb || avgoverfreq || avgovertime || avgoverrpt)
  ft_error('cfg.select ''union'' in combination with averaging across one of the dimensions is not possible');
end

% trim the selection to all inputs, rpt and rpttap are dealt with later
if hasspike,   [selspike,   cfg] = getselection_spike  (cfg, varargin{:}); end
if haspos,     [selpos,     cfg] = getselection_pos    (cfg, varargin{:}, cfg.tolerance, cfg.select); end
if haschan,    [selchan,    cfg] = getselection_chan   (cfg, varargin{:}, cfg.select); end
if haschancmb, [selchancmb, cfg] = getselection_chancmb(cfg, varargin{:}, cfg.select); end
if hasfreq,    [selfreq,    cfg] = getselection_freq   (cfg, varargin{:}, cfg.tolerance, cfg.select); end
if hastime,    [seltime,    cfg] = getselection_time   (cfg, varargin{:}, cfg.tolerance, cfg.select); end

% this is to keep track of all fields that should be retained in the output
keepfield = datfield;

for i=1:numel(varargin)
  
  for j=1:numel(datfield)
    dimtok = tokenize(dimord{j}, '_');
    
    % the rpt selection should only work with a single data argument
    % in case tapers were kept, selrpt~=selrpttap, otherwise selrpt==selrpttap
    [selrpt{i}, dum, rptdim{i}, selrpttap{i}] = getselection_rpt(cfg, varargin{i}, dimord{j});
    
    % check for the presence of each dimension in each datafield
    fieldhasspike   = ismember('spike',   dimtok);
    fieldhaspos     = ismember('pos',     dimtok) || ismember('{pos}', dimtok);
    fieldhaschan    = (ismember('chan',    dimtok) || ismember('{chan}', dimtok)) && isfield(varargin{1}, 'label');
    fieldhaschancmb = ismember('chancmb', dimtok);
    fieldhastime    = ismember('time',    dimtok) && ~hasspike;
    fieldhasfreq    = ismember('freq',    dimtok);
    fieldhasrpt     = ismember('rpt',     dimtok) | ismember('subj', dimtok) | ismember('{rpt}', dimtok);
    fieldhasrpttap  = ismember('rpttap',  dimtok);
    
    % cfg.latency is used to select individual spikes, not to select from a continuously sampled time axis
    
    if fieldhasspike,   varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'spike')),             selspike{i},   false,           'intersect', average); end
    if fieldhaspos,     varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(ismember(dimtok, {'pos', '{pos}'})), selpos{i},     avgoverpos,      cfg.select, average);  end
    if fieldhaschan,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(ismember(dimtok,{'chan' '{chan}'})), selchan{i},    avgoverchan,     cfg.select, average);  end
    if fieldhaschancmb, varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'chancmb')),           selchancmb{i}, avgoverchancmb,  cfg.select, average);  end
    if fieldhastime,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'time')),              seltime{i},    avgovertime,     cfg.select, average);  end
    if fieldhasfreq,    varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, find(strcmp(dimtok,'freq')),              selfreq{i},    avgoverfreq,     cfg.select, average);  end
    if fieldhasrpt,     varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, rptdim{i},                                selrpt{i},     avgoverrpt,      'intersect', average); end
    if fieldhasrpttap,  varargin{i} = makeselection(varargin{i}, datfield{j}, dimtok, rptdim{i},                                selrpttap{i},  avgoverrpt,      'intersect', average); end
    
    % update the fields that should be kept in the structure as a whole
    % and update the dimord for this specific datfield
    keepdim = true(size(dimtok));
    
    if avgoverchan && ~keepchandim
      keepdim(strcmp(dimtok, 'chan')) = false;
      keepfield = setdiff(keepfield, 'label');
    else
      keepfield = [keepfield 'label'];
    end
    
    if avgoverchancmb && ~keepchancmbdim
      keepdim(strcmp(dimtok, 'chancmb')) = false;
      keepfield = setdiff(keepfield, 'labelcmb');
    else
      keepfield = [keepfield 'labelcmb'];
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
      keepdim(strcmp(dimtok, 'dim'))   = false;
      keepfield = setdiff(keepfield, {'pos' '{pos}' 'dim'});
    elseif avgoverpos && keepposdim
      keepfield = setdiff(keepfield, {'dim'}); % this should be removed anyway
    else
      keepfield = [keepfield {'pos' '{pos}' 'dim'}];
    end
    
    if avgoverrpt && ~keeprptdim
      keepdim(strcmp(dimtok, 'rpt'))    = false;
      keepdim(strcmp(dimtok, 'rpttap')) = false;
      keepdim(strcmp(dimtok, 'subj'))   = false;
    end
    
    % update the sampleinfo, if possible, and needed
    if strcmp(datfield{j}, 'sampleinfo') && ~isequal(cfg.latency, 'all')
      if iscell(seltime{i}) && numel(seltime{i})==size(varargin{i}.sampleinfo,1)
        for k = 1:numel(seltime{i})
          varargin{i}.sampleinfo(k,:) = varargin{i}.sampleinfo(k,1) - 1 + seltime{i}{k}([1 end]);
        end
      elseif ~iscell(seltime{i}) && ~isempty(seltime{i}) && ~all(isnan(seltime{i}))
        nrpt       = size(varargin{i}.sampleinfo,1);
        seltime{i} = seltime{i}(:)';
        varargin{i}.sampleinfo = varargin{i}.sampleinfo(:,[1 1]) - 1 + repmat(seltime{i}([1 end]),nrpt,1);
      end
    end
    
    varargin{i}.(datfield{j}) = squeezedim(varargin{i}.(datfield{j}), ~keepdim);
    
  end % for datfield
  
  % also update the fields that describe the dimensions, time/freq/pos have been dealt with as data
  if haschan,    varargin{i} = makeselection_chan   (varargin{i}, selchan{i}, avgoverchan); end % update the label field
  if haschancmb, varargin{i} = makeselection_chancmb(varargin{i}, selchancmb{i}, avgoverchancmb); end % update the labelcmb field
  
end % for varargin

if strcmp(cfg.select, 'union')
  % create the union of the descriptive axes
  if haspos,      varargin = makeunion(varargin, 'pos'); end
  if haschan,     varargin = makeunion(varargin, 'label'); end
  if haschancmb,  varargin = makeunion(varargin, 'labelcmb'); end
  if hastime,     varargin = makeunion(varargin, 'time'); end
  if hasfreq,     varargin = makeunion(varargin, 'freq'); end
end

% remove all fields from the data structure that do not pertain to the selection
sel = strcmp(keepfield, '{pos}'); if any(sel), keepfield(sel) = {'pos'}; end
sel = strcmp(keepfield, 'chan');  if any(sel), keepfield(sel) = {'label'}; end
sel = strcmp(keepfield, 'chancmb');  if any(sel), keepfield(sel) = {'labelcmb'}; end

if avgoverrpt
  % these are invalid after averaging
  keepfield = setdiff(keepfield, {'cumsumcnt' 'cumtapcnt' 'trialinfo' 'sampleinfo'});
end

if avgovertime
  % these are invalid after averaging or making a latency selection
  keepfield = setdiff(keepfield, {'sampleinfo'});
end

for i=1:numel(varargin)
  varargin{i} = keepfields(varargin{i}, [keepfield ignorefields('selectdata')']);
end

% restore the original dimord fields in the data
for i=1:length(orgdim1)
  dimtok = tokenize(orgdim2{i}, '_');
  
  % using a setdiff may result in double occurrences of chan and pos to
  % disappear, so this causes problems as per bug 2962
  % if ~keeprptdim, dimtok = setdiff(dimtok, {'rpt' 'rpttap' 'subj'}); end
  % if ~keepposdim, dimtok = setdiff(dimtok, {'pos' '{pos}'}); end
  % if ~keepchandim, dimtok = setdiff(dimtok, {'chan'}); end
  % if ~keepfreqdim, dimtok = setdiff(dimtok, {'freq'}); end
  % if ~keeptimedim, dimtok = setdiff(dimtok, {'time'}); end
  if ~keeprptdim, dimtok = dimtok(~ismember(dimtok, {'rpt' 'rpttap' 'subj'})); end
  if ~keepposdim, dimtok = dimtok(~ismember(dimtok, {'pos' '{pos}'}));         end
  if ~keepchandim, dimtok = dimtok(~ismember(dimtok, {'chan'})); end
  if ~keepfreqdim, dimtok = dimtok(~ismember(dimtok, {'freq'})); end
  if ~keeptimedim, dimtok = dimtok(~ismember(dimtok, {'time'})); end
  
  dimord = sprintf('%s_', dimtok{:});
  dimord = dimord(1:end-1); % remove the trailing _
  for j=1:length(varargin)
    varargin{j}.(orgdim1{i}) = dimord;
  end
end

% restore the source.avg field, this keeps the output reasonably consistent with the
% old-style source representation of the input
if strcmp(dtype, 'source') && ~isempty(restoreavg)
  for i=1:length(varargin)
    varargin{i}.avg = keepfields(varargin{i}, restoreavg);
    varargin{i}     = removefields(varargin{i}, restoreavg);
  end
end

varargout = varargin;

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance varargout
ft_postamble history varargout
ft_postamble savevar varargout

% the varargout variable can be cleared when written to outputfile
if exist('varargout', 'var') && ft_nargout>numel(varargout)
  % also return the input cfg with the combined selection over all input data structures
  varargout{end+1} = cfg;
end

end % function ft_selectdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection(data, datfield, dimtok, seldim, selindx, avgoverdim, selmode, average)

if numel(seldim) > 1
  for k = 1:numel(seldim)
    data = makeselection(data, datfield, dimtok, seldim(k), selindx, avgoverdim, selmode, average);
  end
  return;
end

if isnumeric(data.(datfield))
  if isrow(data.(datfield)) && seldim==1
    if length(dimtok)==1
      seldim = 2; % switch row and column
    end
  elseif iscolumn(data.(datfield)) && seldim==2
    if length(dimtok)==1
      seldim = 1; % switch row and column
    end
  end
elseif iscell(data.(datfield))
  if isrow(data.(datfield){1}) && seldim==2
    if length(dimtok)==2
      seldim = 3; % switch row and column
    end
  elseif iscolumn(data.(datfield){1}) && seldim==3
    if length(dimtok)==2
      seldim = 2; % switch row and column
    end
  end
end

% an empty selindx means that nothing(!) should be selected and hence everything should be removed, which is different than keeping everything
% the selindx value of NaN indicates that it is not needed to make a selection

switch selmode
  case 'intersect'
    if iscell(selindx)
      % there are multiple selections in multipe vectors, the selection is in the matrices contained within the cell-array
      for j=1:numel(selindx)
        if ~isempty(selindx{j}) && all(isnan(selindx{j}))
          % no selection needs to be made
        else
          data.(datfield){j} = cellmatselect(data.(datfield){j}, seldim-1, selindx{j}, numel(dimtok)==1);
        end
      end
      
    else
      % there is a single selection in a single vector
      if ~isempty(selindx) && all(isnan(selindx))
        % no selection needs to be made
      else
        data.(datfield) = cellmatselect(data.(datfield), seldim, selindx, numel(dimtok)==1);
      end
    end
    
    if avgoverdim
      data.(datfield) = cellmatmean(data.(datfield), seldim, average);
    end
    
  case 'union'
    if ~isempty(selindx) && all(isnan(selindx))
      % no selection needs to be made
    else
      tmp = data.(datfield);
      siz = size(tmp);
      siz(seldim) = numel(selindx);
      data.(datfield) = nan(siz);
      sel = isfinite(selindx);
      switch seldim
        case 1
          data.(datfield)(sel,:,:,:,:,:) = tmp(selindx(sel),:,:,:,:,:);
        case 2
          data.(datfield)(:,sel,:,:,:,:) = tmp(:,selindx(sel),:,:,:,:);
        case 3
          data.(datfield)(:,:,sel,:,:,:) = tmp(:,:,selindx(sel),:,:,:);
        case 4
          data.(datfield)(:,:,:,sel,:,:) = tmp(:,:,:,selindx(sel),:,:);
        case 5
          data.(datfield)(:,:,:,:,sel,:) = tmp(:,:,:,:,selindx(sel),:);
        case 6
          data.(datfield)(:,:,:,:,:,sel) = tmp(:,:,:,:,:,selindx(sel));
        otherwise
          ft_error('unsupported dimension (%d) for making a selection for %s', seldim, datfield);
      end
    end
    
    if avgoverdim
      data.(datfield) = average(data.(datfield), seldim);
    end
end % switch

end % function makeselection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_chan(data, selchan, avgoverchan)
if isempty(selchan)
  %error('no channels were selected');
  data.label = {};
elseif avgoverchan && all(isnan(selchan))
  str = sprintf('%s, ', data.label{:});
  str = str(1:end-2);
  str = sprintf('mean(%s)', str);
  data.label = {str};
elseif avgoverchan && ~any(isnan(selchan))
  str = sprintf('%s, ', data.label{selchan});
  str = str(1:end-2);
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
else
  % this should never happen
  ft_error('cannot figure out how to select channels');
end
end % function makeselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = makeselection_chancmb(data, selchancmb, avgoverchancmb)
if isempty(selchancmb)
  ft_error('no channel combinations were selected');
elseif avgoverchancmb && all(isnan(selchancmb))
  % naming the channel combinations becomes ambiguous, but should not
  % suggest that the mean was computed prior to combining
  str1 = sprintf('%s, ', data.labelcmb{:,1});
  str1 = str1(1:end-2);
  % str1 = sprintf('mean(%s)', str1);
  str2 = sprintf('%s, ', data.labelcmb{:,2});
  str2 = str2(1:end-2);
  % str2 = sprintf('mean(%s)', str2);
  data.label = {str1, str2};
elseif avgoverchancmb && ~any(isnan(selchancmb))
  % naming the channel combinations becomes ambiguous, but should not
  % suggest that the mean was computed prior to combining
  str1 = sprintf('%s, ', data.labelcmb{selchancmb,1});
  str1 = str1(1:end-2);
  % str1 = sprintf('mean(%s)', str1);
  str2 = sprintf('%s, ', data.labelcmb{selchancmb,2});
  str2 = str2(1:end-2);
  % str2 = sprintf('mean(%s)', str2);
  data.label = {str1, str2};
elseif all(isfinite(selchancmb))
  data.labelcmb = data.labelcmb(selchancmb,:);
elseif numel(selchancmb)==1 && any(~isfinite(selchancmb))
  % do nothing
elseif numel(selchancmb)>1  && any(~isfinite(selchancmb))
  tmp = cell(numel(selchancmb),2);
  for k = 1:size(tmp,1)
    if isfinite(selchan(k))
      tmp{k,1} = data.labelcmb{selchan(k),1};
      tmp{k,2} = data.labelcmb{selchan(k),2};
    end
  end
  data.labelcmb = tmp;
else
  % this should never happen
  ft_error('cannot figure out how to select channelcombinations');
end
end % function makeselection_chancmb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chanindx, cfg] = getselection_chan(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

% loop over data once to initialize
chanindx = cell(ndata,1);
label    = cell(1,0);

for k = 1:ndata
  if isfield(varargin{k}, 'grad') && isfield(varargin{k}.grad, 'type')
    % this makes channel selection more robust
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label, varargin{k}.grad.type);
  elseif isfield(varargin{k}, 'elec') && isfield(varargin{k}.elec, 'type')
    % this makes channel selection more robust
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label, varargin{k}.elec.type);
  else
    selchannel = ft_channelselection(cfg.channel, varargin{k}.label);
  end
  label      = union(label, selchannel);
end
label = label(:);   % ensure column array

% this call to match_str ensures that that labels are always in the
% order of the first input argument see bug_2917, but also temporarily keep
% the labels from the other data structures not present in the first one
% (in case selmode = 'union')
[ix, iy] = match_str(varargin{1}.label, label);
label1   = varargin{1}.label(:); % ensure column array
label    = [label1(ix); label(setdiff(1:numel(label),iy))];

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

ok = false(size(indx,1),1);
for k = 1:ndata
  % loop through the columns to preserve the order of the channels, where
  % the order of the input arguments determines the final order
  ix = find(~ok);
  [srt,srtix] = sort(indx(ix,k));
  indx(ix,:)  = indx(ix(srtix),:);
  ok = ok | isfinite(indx(:,k));
end

for k = 1:ndata
  % do a sanity check on double occurrences
  if numel(unique(indx(isfinite(indx(:,k)),k)))<sum(isfinite(indx(:,k)))
    ft_error('the selection of channels across input arguments leads to double occurrences');
  end
  chanindx{k} = indx(:,k);
end

for k = 1:ndata
  if isequal(chanindx{k}, (1:numel(varargin{k}.label))')
    % no actual selection is needed for this data structure
    chanindx{k} = nan;
  end
end

cfg.channel = label;

end % function getselection_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chancmbindx, cfg] = getselection_chancmb(cfg, varargin)

selmode  = varargin{end};
ndata    = numel(varargin)-1;
varargin = varargin(1:ndata);

chancmbindx = cell(ndata,1);

if ~isfield(cfg, 'channelcmb')
  for k=1:ndata
    % the nan return value specifies that no selection was specified
    chancmbindx{k} = nan;
  end
  
else
  
  switch selmode
    case 'intersect'
      for k=1:ndata
        if ~isfield(varargin{k}, 'label')
          cfg.channelcmb = ft_channelcombination(cfg.channelcmb, unique(varargin{k}.labelcmb(:)));
        else
          cfg.channelcmb = ft_channelcombination(cfg.channelcmb, varargin{k}.label);
        end
      end
      
      ncfgcmb = size(cfg.channelcmb,1);
      cfgcmb  = cell(ncfgcmb, 1);
      for i=1:ncfgcmb
        cfgcmb{i} = sprintf('%s&%s', cfg.channelcmb{i,1}, cfg.channelcmb{i,2});
      end
      
      for k=1:ndata
        ndatcmb = size(varargin{k}.labelcmb,1);
        datcmb = cell(ndatcmb, 1);
        for i=1:ndatcmb
          datcmb{i} = sprintf('%s&%s', varargin{k}.labelcmb{i,1}, varargin{k}.labelcmb{i,2});
        end
        
        % return the order according to the (joint) configuration, not according to the (individual) data
        % FIXME this should adhere to the general code guidelines, where
        % the order returned will be according to the first data argument!
        [dum, chancmbindx{k}] = match_str(cfgcmb, datcmb);
      end
      
    case 'union'
      % FIXME this is not yet implemented
      ft_error('union of channel combination is not yet supported');
      
    otherwise
      ft_error('invalid value for cfg.select');
  end % switch
  
  
end

end % function getselection_chancmb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeindx, cfg] = getselection_spike(cfg, varargin)
% possible specifications are
% cfg.latency = string -> 'all'
% cfg.latency = [beg end]
% cfg.trials  = string -> 'all'
% cfg.trials  = vector with indices

ndata    = numel(varargin);
varargin = varargin(1:ndata);

if isequal(cfg.latency, 'all') && isequal(cfg.trials, 'all')
  spikeindx = cell(1,ndata);
  for i=1:ndata
    spikeindx{i} = num2cell(nan(1, length(varargin{i}.time)));
  end
  return
end

trialbeg = varargin{1}.trialtime(:,1);
trialend = varargin{1}.trialtime(:,2);
for i=2:ndata
  trialbeg = cat(1, trialbeg, varargin{1}.trialtime(:,1));
  trialend = cat(1, trialend, varargin{1}.trialtime(:,2));
end

% convert string into a numeric selection
if ischar(cfg.latency)
  switch cfg.latency
    case 'all'
      cfg.latency = [-inf inf];
    case 'maxperiod'
      cfg.latency = [min(trialbeg) max(trialend)];
    case 'minperiod'
      cfg.latency = [max(trialbeg) min(trialend)];
    case 'prestim'
      cfg.latency = [min(trialbeg) 0];
    case 'poststim'
      cfg.latency = [0 max(trialend)];
    otherwise
      ft_error('incorrect specification of cfg.latency');
  end % switch
end

spikeindx = cell(1,ndata);
for i=1:ndata
  nchan = length(varargin{i}.time);
  spikeindx{i} = cell(1,nchan);
  for j=1:nchan
    selbegtime = varargin{i}.time{j}>=cfg.latency(1);
    selendtime = varargin{i}.time{j}<=cfg.latency(2);
    if isequal(cfg.trials, 'all')
      seltrial = true(size(varargin{i}.trial{j}));
    else
      seltrial = ismember(varargin{i}.trial{j}, cfg.trials);
    end
    spikeindx{i}{j} = find(selbegtime & selendtime & seltrial);
  end
end

end % function getselection_spiketime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeindx, cfg] = getselection_time(cfg, varargin)
% possible specifications are
% cfg.latency = value     -> can be 'all'
% cfg.latency = [beg end]

if ft_datatype(varargin{1}, 'spike')
  ft_error('latency selection in spike data is not supported')
end

selmode  = varargin{end};
tol      = varargin{end-1};
ndata    = numel(varargin)-2;
varargin = varargin(1:ndata);

if isequal(cfg.latency, 'all') && iscell(varargin{1}.time)
  % for raw data this means that all trials should be selected as they are
  % for timelock/freq data it is still needed to make the intersection between data arguments
  timeindx = cell(1,ndata);
  for i=1:ndata
    % the nan return value specifies that no selection was specified
    timeindx{i} = num2cell(nan(1, length(varargin{i}.time)));
  end
  return
end

% if there is a single timelock/freq input, there is one time vector
% if there are multiple timelock/freq inputs, there are multiple time vectors
% if there is a single raw input, there are multiple time vectors
% if there are multiple raw inputs, there are multiple time vectors

% collect all time axes in one large cell-array
alltimecell = {};
if iscell(varargin{1}.time)
  for k = 1:ndata
    alltimecell = [alltimecell varargin{k}.time{:}];
  end
else
  for k = 1:ndata
    alltimecell = [alltimecell {varargin{k}.time}];
  end
end

% the nan return value specifies that no selection was specified
timeindx = repmat({nan}, size(alltimecell));

% loop over data once to determine the union of all time axes
alltimevec = zeros(1,0);
for k = 1:length(alltimecell)
  alltimevec = union(alltimevec, round(alltimecell{k}/tol)*tol);
end

indx = nan(numel(alltimevec), numel(alltimecell));
for k = 1:numel(alltimecell)
  [dum, ix, iy] = intersect(alltimevec, round(alltimecell{k}/tol)*tol);
  indx(ix,k) = iy;
end

if iscell(varargin{1}.time) && ischar(cfg.latency)&& ~strcmp(cfg.latency, 'minperiod')
  % if the input data arguments are of type 'raw', temporarily set the
  % selmode to union, otherwise the potentially different length trials
  % will be truncated to the shorted epoch, prior to latency selection.
  selmode = 'union';
elseif ischar(cfg.latency) && strcmp(cfg.latency, 'minperiod')
  % enforce intersect
  selmode = 'intersect';
end
switch selmode
  case 'intersect'
    sel        = sum(isfinite(indx),2)==numel(alltimecell);
    indx       = indx(sel,:);
    alltimevec = alltimevec(sel);
  case 'union'
    % don't do a subselection
  otherwise
    ft_error('invalid value for cfg.select');
end

% Note that cfg.toilim handling has been removed, as it was renamed to cfg.latency

% convert a string selection into a numeric selection
if ischar(cfg.latency)
  switch cfg.latency
    case {'all' 'maxperlen' 'maxperiod'}
      cfg.latency = [min(alltimevec) max(alltimevec)];
    case 'prestim'
      cfg.latency = [min(alltimevec) 0];
    case 'poststim'
      cfg.latency = [0 max(alltimevec)];
    case 'minperiod'
      % the time vector has been pruned above
      cfg.latency = [min(alltimevec) max(alltimevec)];
    otherwise
      ft_error('incorrect specification of cfg.latency');
  end % switch
end

% deal with numeric selection
if isempty(cfg.latency)
  for k = 1:numel(alltimecell)
    % FIXME I do not understand this
    % this signifies that all time bins are deselected and should be removed
    timeindx{k} = [];
  end
  
elseif numel(cfg.latency)==1
  % this single value should be within the time axis of each input data structure
  if numel(alltimevec)>1
    tbin = nearest(alltimevec, cfg.latency, true, true); % determine the numerical tolerance
  else
    tbin = nearest(alltimevec, cfg.latency, true, false); % don't consider tolerance
  end
  cfg.latency = alltimevec(tbin);
  
  for k = 1:ndata
    timeindx{k} = indx(tbin, k);
  end
  
elseif numel(cfg.latency)==2
  % the [min max] range can be specifed with +inf or -inf, but should
  % at least partially overlap with the time axis of the input data
  mintime = min(alltimevec);
  maxtime = max(alltimevec);
  if all(cfg.latency<mintime) || all(cfg.latency>maxtime)
    ft_error('the selected time range falls outside the time axis in the data');
  end
  tbeg = nearest(alltimevec, cfg.latency(1), false, false);
  tend = nearest(alltimevec, cfg.latency(2), false, false);
  cfg.latency = alltimevec([tbeg tend]);
  
  for k = 1:numel(alltimecell)
    timeindx{k} = indx(tbeg:tend, k);
    % if the input data arguments are of type 'raw', the non-finite values
    % need to be removed from the individual cells to ensure correct
    % behavior
    if iscell(varargin{1}.time)
      timeindx{k} = timeindx{k}(isfinite(timeindx{k}));
    end
  end
  
elseif size(cfg.latency,2)==2
  % this may be used for specification of the computation, not for data selection
  
else
  ft_error('incorrect specification of cfg.latency');
end

for k = 1:numel(alltimecell)
  if ~iscell(varargin{1}.time)
    if isequal(timeindx{k}(:)', 1:length(alltimecell{k}))
      % no actual selection is needed for this data structure
      timeindx{k} = nan;
    end
  else
    % if the input data arguments are of type 'raw', they need to be
    % handled differently, because the individual trials can be of
    % different length
  end
end

if iscell(varargin{1}.time)
  % split all time axes again over the different input raw data structures
  dum = cell(1,ndata);
  for k = 1:ndata
    sel = 1:length(varargin{k}.time);
    dum{k} = timeindx(sel); % get the first selection
    timeindx(sel) = []; % remove the first selection
  end
  timeindx = dum;
else
  % no splitting is needed, each input data structure has one selection
end

end % function getselection_time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freqindx, cfg] = getselection_freq(cfg, varargin)
% possible specifications are
% cfg.frequency = value     -> can be 'all'
% cfg.frequency = [beg end]

selmode  = varargin{end};
tol      = varargin{end-1};
ndata    = numel(varargin)-2;
varargin = varargin(1:ndata);

% loop over data once to initialize
freqindx = cell(ndata,1);
freqaxis = zeros(1,0);
for k = 1:ndata
  % the nan return value specifies that no selection was specified
  freqindx{k} = nan;
  
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
  % some of these do not make sense, but are here for consistency with ft_multiplotER
  if ischar(cfg.frequency)
    if strcmp(cfg.frequency, 'all')
      cfg.frequency = [min(freqaxis) max(freqaxis)];
    elseif strcmp(cfg.frequency, 'maxmin')
      cfg.frequency = [min(freqaxis) max(freqaxis)]; % the same as 'all'
    elseif strcmp(cfg.frequency, 'minzero')
      cfg.frequency = [min(freqaxis) 0];
    elseif strcmp(cfg.frequency, 'maxabs')
      cfg.frequency = [-max(abs(freqaxis)) max(abs(freqaxis))];
    elseif strcmp(cfg.frequency, 'zeromax')
      cfg.frequency = [0 max(freqaxis)];
    elseif strcmp(cfg.frequency, 'zeromax')
      cfg.frequency = [0 max(freqaxis)];
    else
      ft_error('incorrect specification of cfg.frequency');
    end
  end
  
  % deal with numeric selection
  if isempty(cfg.frequency)
    for k = 1:ndata
      % FIXME I do not understand this
      % this signifies that all frequency bins are deselected and should be removed
      freqindx{k} = [];
    end
    
  elseif numel(cfg.frequency)==1
    % this single value should be within the frequency axis of each input data structure
    if numel(freqaxis)>1
      fbin = nearest(freqaxis, cfg.frequency, true, true); % determine the numerical tolerance
    else
      fbin = nearest(freqaxis, cfg.frequency, true, false); % don't consider tolerance
    end
    cfg.frequency = freqaxis(fbin);
    
    for k = 1:ndata
      freqindx{k} = indx(fbin,k);
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
      freqindx{k} = indx(fbeg:fend,k);
    end
    
  elseif size(cfg.frequency,2)==2
    % this may be used for specification of the computation, not for data selection
    
  else
    ft_error('incorrect specification of cfg.frequency');
  end
end % if cfg.frequency

for k = 1:ndata
  if isequal(freqindx{k}, 1:length(varargin{k}.freq))
    % the cfg was updated, but no selection is needed for the data
    freqindx{k} = nan;
  end
end

end % function getselection_freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rptindx, cfg, rptdim, rpttapindx] = getselection_rpt(cfg, varargin)
% this should deal with cfg.trials

dimord   = varargin{end};
ndata    = numel(varargin)-1;
data     = varargin{1:ndata}; % this syntax ensures that it will only work on a single data input

dimtok = tokenize(dimord, '_');
rptdim = find(strcmp(dimtok, '{rpt}') | strcmp(dimtok, 'rpt') | strcmp(dimtok, 'rpttap') | strcmp(dimtok, 'subj'));

if isequal(cfg.trials, 'all')
  rptindx    = nan; % the nan return value specifies that no selection was specified
  rpttapindx = nan; % the nan return value specifies that no selection was specified
  
elseif isempty(rptdim)
  % FIXME should [] not mean that none of the trials is to be selected?
  rptindx    = nan; % the nan return value specifies that no selection was specified
  rpttapindx = nan; % the nan return value specifies that no selection was specified
  
else
  rptindx = ft_getopt(cfg, 'trials');
  
  if islogical(rptindx)
    % convert from booleans to indices
    rptindx = find(rptindx);
  end
  
  rptindx = unique(sort(rptindx));
  
  if strcmp(dimtok{rptdim}, 'rpttap') && isfield(data, 'cumtapcnt')
    % there are tapers in the data
    
    % determine for each taper to which trial it belongs
    nrpt = size(data.cumtapcnt, 1);
    taper = zeros(nrpt, 1);
    sumtapcnt = cumsum([0; data.cumtapcnt(:)]);
    begtapcnt = sumtapcnt(1:end-1)+1;
    endtapcnt = sumtapcnt(2:end);
    for i=1:nrpt
      taper(begtapcnt(i):endtapcnt(i)) = i;
    end
    rpttapindx = find(ismember(taper, rptindx));
    
  else
    % there are no tapers in the data
    rpttapindx = rptindx;
  end
end

end % function getselection_rpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posindx, cfg] = getselection_pos(cfg, varargin)
% possible specifications are <none>

ndata   = numel(varargin)-2;
tol     = varargin{end-1}; % FIXME this is still ignored
selmode = varargin{end};   % FIXME this is still ignored
data    = varargin(1:ndata);

for i=1:ndata
  if ~isequal(varargin{i}.pos, varargin{1}.pos)
    % FIXME it would be possible here to make a selection based on intersect or union
    ft_error('not yet implemented');
  end
end

if strcmp(cfg.select, 'union')
  % FIXME it would be possible here to make a selection based on intersect or union
  ft_error('not yet implemented');
end

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
if isvector(x)
  % there is no harm to keep it as it is
elseif istable(x)
  % there is no harm to keep it as it is
else
  x = reshape(x, [siz(~dim) 1]);
end
end % function squeezedim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = makeunion(x, field)
old = cellfun(@getfield, x, repmat({field}, size(x)), 'uniformoutput', false);
if iscell(old{1})
  % empty is indicated to represent missing value for a cell-array (label, labelcmb)
  new = old{1};
  for i=2:length(old)
    sel = ~cellfun(@isempty, old{i});
    new(sel) = old{i}(sel);
  end
else
  % nan is indicated to represent missing value for a numeric array (time, freq, pos)
  new = old{1};
  for i=2:length(old)
    sel = ~isnan(old{i});
    new(sel) = old{i}(sel);
  end
end
x = cellfun(@setfield, x, repmat({field}, size(x)), repmat({new}, size(x)), 'uniformoutput', false);
end % function makeunion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to make a selextion in data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cellmatselect(x, seldim, selindx, maybevector)
if nargin<4
  % some fields are a vector with an unspecified singleton dimension, these can be transposed
  % if the singleton dimension represents something explicit, they should not be transposed
  % they might for example represent a single trial, or a single channel
  maybevector = true;
end
if iscell(x)
  if seldim==1
    x = x(selindx);
  else
    for i=1:numel(x)
      if isempty(x{i})
        continue
      end
      switch seldim
        case 2
          if maybevector && isvector(x{i})
            % sometimes the data is 1xN, whereas the dimord describes only the first dimension
            % in this case a row and column vector can be interpreted as equivalent
            x{i} = x{i}(selindx);
          elseif istable(x)
            % multidimensional indexing is not supported
            x{i} = x{i}(selindx,:);
          else
            x{i} = x{i}(selindx,:,:,:,:);
          end
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
      if maybevector && isvector(x)
        % sometimes the data is 1xN, whereas the dimord describes only the first dimension
        % in this case a row and column vector can be interpreted as equivalent
        x = x(selindx);
      elseif istable(x)
        % multidimensional indexing is not supported
        x = x(selindx,:);
      else
        x = x(selindx,:,:,:,:,:);
      end
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
function x = cellmatmean(x, seldim, average)
if iscell(x)
  if seldim==1
    for i=2:numel(x)
      x{1} = x{1} + x{i};
    end
    x = {x{1}/numel(x)};
  else
    for i=1:numel(x)
      x{i} = average(x{i}, seldim-1);
    end % for
  end
else
  x = average(x, seldim);
end
end % function cellmatmean
