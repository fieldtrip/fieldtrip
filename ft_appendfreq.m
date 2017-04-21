function [freq] = ft_appendfreq(cfg, varargin)

% FT_APPENDFREQ concatenates multiple frequency or time-frequency data structures
% that have been processed separately. If the input data structures contain different
% channels, it will be concatenated along the channel direction. If the channels are
% identical in the input data structures, the data will be concatenated along the
% repetition dimension.
%
% Use as
%  combined = ft_appendfreq(cfg, freq1, freq2, ...)
%
% The configuration should contain
%   cfg.parameter  = string, the name of the field to concatenate
%
% The configuration can optionally contain
%   cfg.appenddim  = string, the dimension to concatenate over (default is automatic)
%   cfg.tolerance  = scalar, tolerance to determine how different the frequency and/or
%                    time axes are allowed to still be considered compatible (default = 1e-5)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a
% *.mat file on disk and/or the output data will be written to a *.mat file.
% These mat files should contain only a single variable, corresponding with
% the input/output structure.
%
% See also FT_FREQANALYSIS, FT_DATATYPE_FREQ, FT_APPENDDATA, FT_APPENDTIMELOCK,
% FT_APPENDSENS

% Copyright (C) 2011-2017, Robert Oostenveld
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

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  % FIXME: what about freq+comp?
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel', 'all');
cfg.parameter  = ft_getopt(cfg, 'parameter', []);
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5);
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');

hastime = isfield(varargin{1}, 'time');
hasfreq = isfield(varargin{1}, 'freq');

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if hastime && hasfreq
    if checkchan(varargin{:}, 'identical') && checkfreq(varargin{:}, 'identical', cfg.tolerance) && checktime(varargin{:}, 'identical', cfg.tolerance)
      cfg.appenddim = 'rpt';
    elseif checkfreq(varargin{:}, 'unique', cfg.tolerance) && checktime(varargin{:}, 'identical', cfg.tolerance)
      cfg.appenddim = 'freq';
    elseif checktime(varargin{:}, 'unique', cfg.tolerance) && checkfreq(varargin{:}, 'identical', cfg.tolerance)
      cfg.appenddim = 'time';
    elseif checkchan(varargin{:}, 'unique')
      cfg.appenddim = 'chan';
    else
      error('cfg.appenddim should be specified');
    end
  else
    if checkchan(varargin{:}, 'identical') && checkfreq(varargin{:}, 'identical', cfg.tolerance)
      cfg.appenddim = 'rpt';
    elseif checkfreq(varargin{:}, 'unique', cfg.tolerance)
      cfg.appenddim = 'freq';
    elseif checkchan(varargin{:}, 'unique')
      cfg.appenddim = 'chan';
    else
      error('cfg.appenddim should be specified');
    end
  end
end

fprintf('concatenating over the "%s" dimension\n', cfg.appenddim);

if isempty(cfg.parameter)
  fn = fieldnames(varargin{1});
  for i=2:numel(varargin)
    fn = intersect(fn, fieldnames(varargin{i}));
  end
  fn = setdiff(fn, ignorefields('appendtimelock'));
  assert(~isempty(fn), 'cfg.parameter shoudl be specified');
else
  fn = {cfg.parameter};
end

switch cfg.appenddim
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'chan'
    assert(checkchan(varargin{:}, 'unique'));
    % remember the original channel labels in each input
    oldlabel = cell(size(varargin));
    for i=1:numel(varargin)
      oldlabel{i} =  varargin{i}.label;
    end
    
    % determine the union of all input data
    tmpcfg = [];
    tmpcfg.select = 'union';
    tmpcfg.tolerance = cfg.tolerance;
    tmpcfg.channel = cfg.channel;
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the union of all input data
    freq = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord'});
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_freq'
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            freq.(fn{i})(chansel,:) = varargin{j}.(fn{i})(chansel,:);
          end
          
        case {'rpt_chan_freq' 'subj_chan_freq'}
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            freq.(fn{i})(:,chansel,:) = varargin{j}.(fn{i})(:,chansel,:);
          end
          
        otherwise
          % do not concatenate this field
      end % switch
    end % for fn
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'freq'
    assert(checkfreq(varargin{:}, 'unique', cfg.tolerance));
    
    % remember the original time axes in each input
    oldfreq = cell(size(varargin));
    for i=1:numel(varargin)
      oldfreq{i} =  varargin{i}.freq;
    end
    
    % determine the union of all input data
    tmpcfg = [];
    tmpcfg.select = 'union';
    tmpcfg.tolerance = cfg.tolerance;
    tmpcfg.channel = cfg.channel;
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the union of all input data
    freq = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord'});
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_freq'
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            freq.(fn{i})(:,freqsel) = varargin{j}.(fn{i})(:,freqsel);
          end
          
        case 'chan_freq_time'
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            freq.(fn{i})(:,freqsel,:) = varargin{j}.(fn{i})(:,freqsel,:);
          end
          
        case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            freq.(fn{i})(:,:,freqsel) = varargin{j}.(fn{i})(:,:,freqsel);
          end
          
        case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            freq.(fn{i})(:,:,freqsel,:) = varargin{j}.(fn{i})(:,:,freqsel,:);
          end
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for fn
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'time'
    assert(checktime(varargin{:}, 'unique', cfg.tolerance));
    
    % remember the original time axes in each input
    oldtime = cell(size(varargin));
    for i=1:numel(varargin)
      oldtime{i} =  varargin{i}.time;
    end
    
    % determine the union of all input data
    tmpcfg = [];
    tmpcfg.select = 'union';
    tmpcfg.tolerance = cfg.tolerance;
    tmpcfg.channel = cfg.channel;
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the union of all input data
    freq = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord'});
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_freq'
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            freq.(fn{i})(:,timsel) = varargin{j}.(fn{i})(:,timsel);
          end
          
        case 'chan_freq_time'
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            freq.(fn{i})(:,:,timsel) = varargin{j}.(fn{i})(:,:,timsel);
          end
          
        case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            freq.(fn{i})(:,:,timsel) = varargin{j}.(fn{i})(:,:,timsel);
          end
          
        case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
          freq.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            freq.(fn{i})(:,:,:,timsel) = varargin{j}.(fn{i})(:,:,:,timsel);
          end
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for fn
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'rpt'
    % determine the intersection of all input data
    tmpcfg = [];
    tmpcfg.select = 'intersect';
    tmpcfg.tolerance = cfg.tolerance;
    tmpcfg.channel = cfg.channel;
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the intersection of all input data
    freq = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord'});
    assert(numel(freq.label)>0);
    assert(numel(freq.time)>0);
    assert(numel(freq.freq)>0);
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_freq'
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = reshape(varargin{j}.(fn{i}), 1, dimsiz(1), dimsiz(2));
          end
          freq.(fn{i}) = cat(1, dat{:});
          
        case 'chan_freq_time'
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = reshape(varargin{j}.(fn{i}), 1, dimsiz(1), dimsiz(2), dimsiz(3));
          end
          freq.(fn{i}) = cat(1, dat{:});
          
        case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = varargin{j}.(fn{i});
          end
          freq.(fn{i}) = cat(1, dat{:});
          
        case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = varargin{j}.(fn{i});
          end
          freq.(fn{i}) = cat(1, dat{:});
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for fn
    
  otherwise
    error('unsupported cfg.appenddim');
end

if isfield(freq, 'dimord')
  dimtok = tokenize(freq.dimord);
  if strcmp(cfg.appenddim, 'rpt') && ~any(strcmp(dimtok{1}, {'rpt', 'rpttap', 'subj'}))
    freq.dimord = ['rpt_' freq.dimord];
  end
end

% % create the output structure from scratch
% freq   = [];
%
% tol    = cfg.tolerance;
% dimtok = tokenize(dimord{1}, '_');
%
% if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
%   % only allow to append across observations if these are present in the data
%   if any(strcmp(dimtok, 'rpt'))
%     cfg.appenddim = 'rpt';
%   elseif any(strcmp(dimtok, 'rpttap'))
%     cfg.appenddim = 'rpttap';
%   elseif any(strcmp(dimtok, 'subj'))
%     cfg.appenddim = 'subj';
%   else
%     % we need to check whether the other dimensions are the same.
%     % if not, consider some tolerance.
%     boolval1 = checkchan(varargin{:}, 'identical');
%     boolval2 = checkfreq(varargin{:}, 'identical', tol);
%
%     if isfield(varargin{1}, 'time')
%       boolval3 = checktime(varargin{:}, 'identical', tol);
%       if boolval1 && boolval2 && boolval3
%         % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
%         cfg.appenddim = 'rpt';
%       elseif ~boolval1 && boolval2 && boolval3
%         cfg.appenddim = 'chan';
%       elseif boolval1 && ~boolval2 && boolval3
%         cfg.appenddim = 'freq';
%       elseif boolval1 && boolval2 && ~boolval3
%         cfg.appenddim = 'time';
%       else
%         error('the input datasets have multiple non-identical dimensions, this function only appends one dimension at a time');
%       end
%     else
%       if boolval1 && boolval2
%         % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
%         cfg.appenddim = 'rpt';
%       elseif ~boolval1 && boolval2
%         cfg.appenddim = 'chan';
%       elseif boolval1 && ~boolval2
%         cfg.appenddim = 'freq';
%       end
%     end
%
%   end % determine the dimension for appending
% end % isempty
%
% switch cfg.appenddim
%   case {'rpt' 'rpttap' 'subj'}
%     catdim = find(ismember(dimtok, {'rpt' 'rpttap' 'subj'}));
%     if numel(catdim)==0
%       catdim = 0;
%     elseif numel(catdim)==1
%       % this is OK
%     elseif numel(catdim)>1
%       error('ambiguous dimord for concatenation');
%     end
%
%     % if any of these are present, concatenate
%     % if not prepend the dimord with rpt (and thus shift the dimensions)
%
%     % here we need to check whether the other dimensions are the same. if
%     % not, consider some tolerance.
%     boolval1 = checkchan(varargin{:}, 'identical');
%     boolval2 = checkfreq(varargin{:}, 'identical', tol);
%     if isfield(varargin{1}, 'time')
%       boolval3 = checktime(varargin{:}, 'identical', tol);
%     else
%       boolval3 = true;
%     end
%
%     if any([boolval2 boolval3]==false)
%       error('appending across observations is not possible, because the frequency and/or temporal dimensions are incompatible');
%     end
%
%     % select and reorder the channels that are in every dataset
%     tmpcfg           = [];
%     tmpcfg.channel   = cfg.channel;
%     tmpcfg.tolerance = cfg.tolerance;
%     [varargin{:}]    = ft_selectdata(tmpcfg, varargin{:});
%     for i=1:Ndata
%       [cfg_rolledback, varargin{i}] = rollback_provenance(cfg, varargin{i});
%     end
%     cfg = cfg_rolledback;
%
%     % update the dimord
%     if catdim==0
%       freq.dimord = ['rpt_',varargin{1}.dimord];
%       % FIXME append dof
%     else
%       freq.dimord = varargin{1}.dimord;
%       % FIXME append dof
%       % before append cumtapcnt cumsumcnt trialinfo check if there's a
%       % subfield in each dataset. Append fields of different dataset might
%       % lead in empty and/or non-existing fields in a particular dataset
%       hascumsumcnt = [];
%       hascumtapcnt = [];
%       hastrialinfo = [];
%       for i=1:Ndata
%         if isfield(varargin{i},'cumsumcnt')
%           hascumsumcnt(end+1) = 1;
%         else
%           hascumsumcnt(end+1) = 0;
%         end
%         if isfield(varargin{i},'cumtapcnt')
%           hascumtapcnt(end+1) = 1;
%         else
%           hascumtapcnt(end+1) = 0;
%         end
%         if isfield(varargin{i},'trialinfo')
%           hastrialinfo(end+1) = 1;
%         else
%           hastrialinfo(end+1) = 0;
%         end
%       end
%
%       % screen concatenable fields
%       if ~checkfreq(varargin{:}, 'identical', tol)
%         error('the freq fields of the input data structures are not equal');
%       else
%         freq.freq=varargin{1}.freq;
%       end
%       if ~sum(hascumsumcnt)==0 && ~(sum(hascumsumcnt)==Ndata)
%         error('the cumsumcnt fields of the input data structures are not equal');
%       else
%         iscumsumcnt=unique(hascumsumcnt);
%       end
%       if ~sum(hascumtapcnt)==0 && ~(sum(hascumtapcnt)==Ndata)
%         error('the cumtapcnt fields of the input data structures are not equal');
%       else
%         iscumtapcnt=unique(hascumtapcnt);
%       end
%       if ~sum(hastrialinfo)==0 && ~(sum(hastrialinfo)==Ndata)
%         error('the trialinfo fields of the input data structures are not equal');
%       else
%         istrialinfo=unique(hastrialinfo);
%       end
%
%       % concatenating fields
%       for i=1:Ndata
%         if iscumsumcnt
%           cumsumcnt{i}=varargin{i}.cumsumcnt;
%         end
%         if iscumtapcnt
%           cumtapcnt{i}=varargin{i}.cumtapcnt;
%         end
%         if istrialinfo
%           trialinfo{i}=varargin{i}.trialinfo;
%         end
%       end
%
%       % fill in the rest of the descriptive fields
%       if iscumsumcnt
%         freq.cumsumcnt = cat(catdim,cumsumcnt{:});
%         clear cumsumcnt;
%       end
%       if iscumtapcnt
%         freq.cumtapcnt = cat(catdim,cumtapcnt{:});
%         clear cumtapcnt;
%       end
%       if istrialinfo
%         freq.trialinfo = cat(catdim,trialinfo{:});
%         clear trialinfo;
%       end
%     end
%
%     freq.label = varargin{1}.label;
%     freq.freq  = varargin{1}.freq;
%     if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
%
%   case 'chan'
%     catdim = find(strcmp('chan', dimtok));
%     if isempty(catdim)
%       % try chancmb
%       catdim = find(strcmp('chancmb', dimtok));
%     elseif numel(catdim)>1
%       error('ambiguous dimord for concatenation');
%     end
%
%     % check whether all channels are unique and throw an error if not
%     [boolval, list] = checkchan(varargin{:}, 'unique');
%     if ~boolval
%       error('the input data structures have non-unique channels, concatenation across channel is not possible');
%     end
%
%     if isfield(varargin{1}, 'time')
%       if ~checktime(varargin{:}, 'identical', tol)
%         error('the input data structures have non-identical time bins, concatenation across channels not possible');
%       end
%     end
%
%     if ~checkfreq(varargin{:}, 'identical', tol)
%       error('the input data structures have non-identical frequency bins, concatenation across channels not possible');
%     end
%
%     % update the channel description
%     freq.label = list;
%
%     % fill in the rest of the descriptive fields
%     freq.freq   = varargin{1}.freq;
%     if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
%     freq.dimord = varargin{1}.dimord;
%
%   case 'freq'
%     catdim = find(strcmp('freq', dimtok));
%
%     % check whether all frequencies are unique and throw an error if not
%     [boolval, list] = checkfreq(varargin{:}, 'unique', tol);
%     if ~boolval
%       error('the input data structures have non-unique frequency bins, concatenation across frequency is not possible');
%     end
%
%     if ~checkchan(varargin{:}, 'identical')
%       error('the input data structures have non-identical channels, concatenation across frequency not possible');
%     end
%
%     if isfield(varargin{1}, 'time')
%       if ~checktime(varargin{:}, 'identical', tol)
%         error('the input data structures have non-identical time bins, concatenation across channels not possible');
%       end
%     end
%
%     % update the frequency description
%     freq.freq = list(:)';
%
%     % fill in the rest of the descriptive fields
%     freq.label  = varargin{1}.label;
%     freq.dimord = varargin{1}.dimord;
%     if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
%
%   case 'time'
%     catdim = find(strcmp('time', dimtok));
%
%     % check whether all time points are unique and throw an error if not
%     [boolval, list] = checktime(varargin{:}, 'unique', tol);
%     if ~boolval
%       error('the input data structures have non-unique time bins, concatenation across time is not possible');
%     end
%
%     if ~checkchan(varargin{:}, 'identical')
%       error('the input data structures have non-identical channels, concatenation across time not possible');
%     end
%     if ~checkfreq(varargin{:}, 'identical', tol)
%       error('the input data structures have non-identical frequency bins, concatenation across time not possible');
%     end
%
%     % update the time description
%     freq.time = list(:)';
%
%     % fill in the rest of the descriptive fields
%     freq.label  = varargin{1}.label;
%     freq.freq   = varargin{1}.freq;
%     freq.dimord = varargin{1}.dimord;
%
%   otherwise
%     error('it is not allowed to concatenate across dimension %s',cfg.appenddim);
% end
%
% param = cfg.parameter;
% if ~iscell(param), param = {param}; end
%
% % are we appending along the channel dimension?
% catchan = strcmp(cfg.appenddim, 'chan');
% chandim = find(strcmp('chan', dimtok));
%
% % concatenate the numeric data
% for k = 1:numel(param)
%   tmp = cell(size(varargin));
%
%   % get the numeric data 'param{k}' if present
%   for m = 1:Ndata
%
%     if ~isfield(varargin{m}, param{k})
%       error('parameter %s is not present in all data sets', param{k});
%     end
%     tmp{m} = varargin{m}.(param{k});
%
%     % if we are not appending along the channel dimension, make sure we
%     % reorder the channel dimension across the different data sets. At this
%     % point we can be sure that all data sets have identical channels.
%     if ~catchan && m > 1
%       [a,b] = match_str(varargin{1}.label, varargin{m}.label);
%       if ~all(a==b)
%         tmp{m} = reorderdim(tmp{m}, chandim, b);
%       end
%     end
%   end
%
%   if catdim==0
%     ndim    = length(size(tmp{1}));
%     freq.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
%   else
%     freq.(param{k}) = cat(catdim,tmp{:});
%   end
% end % for k = 1:numel(param)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following section is shared with ft_appenddata, ft_appendtimelock and ft_appendfreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hasgrad = false;
haselec = false;
hasopto = false;
for i=1:length(varargin)
  hasgrad = hasgrad || isfield(varargin{i}, 'grad');
  haselec = haselec || isfield(varargin{i}, 'elec');
  hasopto = hasopto || isfield(varargin{i}, 'opto');
end
if  hasgrad || haselec || hasopto
  % gather the sensor definitions from all inputs
  grad = cell(size(varargin));
  elec = cell(size(varargin));
  opto = cell(size(varargin));
  for j=1:Ndata
    if isfield(varargin{j}, 'elec')
      elec{j} = varargin{j}.elec;
    end
    if isfield(varargin{j}, 'grad')
      grad{j} = varargin{j}.grad;
    end
    if isfield(varargin{j}, 'opto')
      opto{j} = varargin{j}.opto;
    end
  end
  % see test_pull393.m for a description of the expected behavior
  if strcmp(cfg.appendsens, 'yes')
    fprintf('concatenating sensor information across input arguments\n');
    % append the sensor descriptions, skip the empty ones
    if hasgrad, freq.grad = ft_appendsens([], grad{~cellfun(@isempty, grad)}); end
    if haselec, freq.elec = ft_appendsens([], elec{~cellfun(@isempty, elec)}); end
    if hasopto, freq.opto = ft_appendsens([], opto{~cellfun(@isempty, opto)}); end
  else
    % discard sensor information when it is inconsistent across the input arguments
    removegrad = any(cellfun(@isempty, grad));
    removeelec = any(cellfun(@isempty, elec));
    removeopto = any(cellfun(@isempty, opto));
    for j=2:length(varargin)
      removegrad = removegrad || ~isequaln(grad{j}, grad{1});
      removeelec = removeelec || ~isequaln(elec{j}, elec{1});
      removeopto = removeopto || ~isequaln(opto{j}, opto{1});
    end
    if hasgrad && ~removegrad, freq.grad = grad{1}; end
    if haselec && ~removeelec, freq.elec = elec{1}; end
    if hasopto && ~removeopto, freq.opto = opto{1}; end
  end
end


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance freq
ft_postamble history    freq
ft_postamble savevar    freq
