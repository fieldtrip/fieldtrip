function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data structures that
% have been processed seperately. If the input data structures contain different
% channels, it will be concatenated along the channel direction. If the channels are
% identical in the input data structures, the data will be concatenated along the
% repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% The configuration can optionally contain
%   cfg.appenddim  = string, the dimension to concatenate over which to append,
%                    this can be 'chan' and 'rpt' (default is automatic)
%   cfg.tolerance  = scalar, tolerance to determine how different the time axes
%                    are allowed to still be considered compatible (default = 1e-5)
%
% See also FT_TIMELOCKANALYSIS, FT_DATATYPE_TIMELOCK, FT_APPENDDATA, FT_APPENDFREQ,
% FT_APPENDSENS, FT_APPENDSOURCE

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
  % FIXME: what about timelock+comp?
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel', 'all');
cfg.parameter  = ft_getopt(cfg, 'parameter', []);
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5);
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if checkchan(varargin{:}, 'identical') && checktime(varargin{:}, 'identical', cfg.tolerance)
    cfg.appenddim = 'rpt';
  elseif checktime(varargin{:}, 'identical', cfg.tolerance) && ~checkchan(varargin{:}, 'unique')
    cfg.appenddim = 'rpt';
  elseif checkchan(varargin{:}, 'unique')
    cfg.appenddim = 'chan';
  elseif checktime(varargin{:}, 'unique', cfg.tolerance)
    cfg.appenddim = 'time';
  else
    error('cfg.appenddim should be specified');
  end
end

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
    timelock = keepfields(varargin{1}, {'label', 'time'});
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_time'
          timelock.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            timelock.(fn{i})(chansel,:) = varargin{j}.(fn{i})(chansel,:);
          end
          
        case {'rpt_chan_time' 'subj_chan_time'}
          timelock.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            timelock.(fn{i})(:,chansel,:) = varargin{j}.(fn{i})(:,chansel,:);
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
    timelock = keepfields(varargin{1}, {'label', 'time'});
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_time'
          timelock.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            timelock.(fn{i})(:,timsel) = varargin{j}.(fn{i})(:,timsel);
          end
          
        case {'rpt_chan_time' 'subj_chan_time'}
          timelock.(fn{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timsel = match_val(varargin{j}.time, oldtime{j});
            timelock.(fn{i})(:,:,timsel) = varargin{j}.(fn{i})(:,:,timsel);
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
    timelock = keepfields(varargin{1}, {'label', 'time'});
    assert(numel(timelock.label)>0);
    assert(numel(timelock.time)>0);
    
    for i=1:numel(fn)
      dimsiz = getdimsiz(varargin{1}, fn{i});
      switch getdimord(varargin{1}, fn{i})
        case 'chan_time'
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = reshape(varargin{j}.(fn{i}), 1, dimsiz(1), dimsiz(2));
          end
          timelock.(fn{i}) = cat(1, dat{:});
          
        case {'rpt_chan_time' 'subj_chan_time'}
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = varargin{j}.(fn{i});
          end
          timelock.(fn{i}) = cat(1, dat{:});
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for fn
    
  otherwise
    error('unsupported cfg.appenddim');
end

if isfield(timelock, 'dimord')
  dimtok = tokenize(timelock.dimord);
  if strcmp(cfg.appenddim, 'rpt') && ~any(strcmp(dimtok{1}, {'rpt', 'subj'}))
    timelock.dimord = ['rpt_' timelock.dimord];
  end
end

% % ensure consistent input data
% if ~checktime(varargin{:}, 'identical', cfg.tolerance)
%   error('this function requires identical time axes for all input structures');
% end
%
% % do a basic check to see whether the dimords match
% dimord = cell(size(varargin));
% for i=1:Ndata
%   dimord{i} = varargin{i}.dimord;
% end
% dimordmatch = all(strcmp(dimord{1}, dimord));
%
% if ~dimordmatch
%   error('the dimords of the input data structures are not equal');
% end
%
% % determine over which dimension to append if not user-specified
% if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
%   [boolval, list] = checkchan(varargin{:}, 'unique');
%   if boolval
%     cfg.appenddim='chan';
%   else
%     cfg.appenddim='rpt';
%   end
% end % isempty
%
% % start with the initial output structure
% timelock        = [];
% timelock.time   = varargin{1}.time;
% ntime           = length(timelock.time);
%
% switch cfg.appenddim
%   case {'rpt' 'rpttap' 'subj'}
%     % select the channels that are common to every dataset
%     tmpcfg           = [];
%     tmpcfg.channel   = cfg.channel;
%     tmpcfg.tolerance = cfg.tolerance;
%     [varargin{:}]    = ft_selectdata(tmpcfg, varargin{:});
%     for i=1:Ndata
%       [cfg_rolledback, varargin{i}] = rollback_provenance(cfg, varargin{i});
%     end
%     cfg = cfg_rolledback;
%
%     timelock.label = varargin{1}.label;
%     nchan          = numel(timelock.label);
%     if nchan<1
%       error('No channels in common');
%     end
%
%     hascov = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==3;
%     if isfield(varargin{1}, 'trial')
%       % these don't make sense when concatenating the avg
%       hastrialinfo  = isfield(varargin{1}, 'trialinfo');
%       hassampleinfo = isfield(varargin{1}, 'sampleinfo');
%
%       ntrial = zeros(size(varargin));
%       for i=1:length(varargin)
%         ntrial(i) = size(varargin{i}.trial, 1);
%       end
%       trialsel = cumsum([1 ntrial]);
%
%       timelock.trial = zeros(sum(ntrial), nchan, ntime);
%       if hastrialinfo,  timelock.trialinfo = zeros(sum(ntrial), size(varargin{1}.trialinfo,2)); end
%       if hassampleinfo, timelock.sampleinfo = zeros(sum(ntrial), size(varargin{1}.sampleinfo,2)); end
%       if hascov, timelock.cov = zeros(sum(ntrial), nchan, nchan); end
%
%       for i=1:length(varargin)
%         % copy the desired data into the output structure
%         begtrial = trialsel(i);
%         endtrial = trialsel(i+1)-1;
%         chansel = match_str(cfg.channel, varargin{i}.label);
%         timelock.trial(begtrial:endtrial,:,:) = varargin{i}.trial(:,chansel,:);
%         if hastrialinfo,  timelock.trialinfo(begtrial:endtrial,:)   = varargin{i}.trialinfo(:,:); end
%         if hassampleinfo, timelock.sampleinfo(begtrial:endtrial,:)  = varargin{i}.sampleinfo(:,:); end
%         if hascov,        timelock.cov(begtrial:endtrial,:,:)       = varargin{i}.cov(:,chansel,chansel); end
%       end % for varargin
%       timelock.avg = permute(mean(timelock.trial,1),[2 3 1]);
%       timelock.var = permute(var(timelock.trial,0,1),[2 3 1]);
%
%     elseif isfield(varargin{1}, 'avg')
%
%       ntrial = numel(varargin);
%       timelock.trial = zeros(ntrial, nchan, ntime);
%       if hascov, timelock.cov = zeros(sum(ntrial),nchan,nchan); end
%
%       for i=1:length(varargin)
%         % copy the desired data into the output structure
%         chansel = match_str(cfg.channel, varargin{i}.label);
%         timelock.trial(i,:,:) = varargin{i}.avg(chansel,:);
%         if hascov, timelock.cov(i,:,:) = varargin{i}.cov(chansel,chansel); end
%       end % for varargin
%     end
%     timelock.dimord = 'rpt_chan_time';
%
%   case 'chan'
%     for i=2:length(varargin) % check if any channels in common
%       if ~isempty(ft_channelselection(varargin{1}.label, varargin{i}.label))
%         error('Cannot concatenate over channels, since some channels are in common');
%       end
%     end
%     % no channels in common; append over channels
%     nchan = zeros(size(varargin));
%     for i=1:length(varargin)
%       nchan(i) = length(varargin{i}.label);
%     end
%     chansel = cumsum([1 nchan]);
%     timelock.label = cell(0,1);
%
%     hascov        = isfield(varargin{1}, 'cov') && numel(size(varargin{1}.cov))==3;
%     if hascov, warning('Concatenating over channels does not allow for keeping covariance information'); end
%
%     if isfield(varargin{1}, 'trial')
%       timelock.dimord = 'rpt_chan_time';
%
%       % these don't make sense when concatenating the avg
%       hastrialinfo  = isfield(varargin{1}, 'trialinfo');
%       hassampleinfo = isfield(varargin{1}, 'sampleinfo');
%
%       ntrial = size(varargin{1}.trial,1);
%
%       timelock.trial = zeros(ntrial, sum(nchan), ntime);
%       if hastrialinfo,  timelock.trialinfo = varargin{1}.trialinfo; end
%       if hassampleinfo, timelock.sampleinfo = varargin{1}.sampleinfo; end
%
%       for i=1:length(varargin)
%         % copy the desired data into the output structure
%         begchan = chansel(i);
%         endchan = chansel(i+1)-1;
%         timelock.trial(:,begchan:endchan,:) = varargin{i}.trial;
%         timelock.label = [timelock.label; varargin{i}.label(:)];
%       end % for varargin
%       timelock.avg = permute(mean(timelock.trial,1),[2 3 1]);
%       timelock.var = permute(var(timelock.trial,0,1),[2 3 1]);
%       timelock.dimord = 'rpt_chan_time';
%
%     else
%       timelock.avg = zeros(sum(nchan), ntime);
%
%       for i=1:length(varargin)
%         % copy the desired data into the output structure
%         begchan = chansel(i);
%         endchan = chansel(i+1)-1;
%         timelock.avg(begchan:endchan,:) = varargin{i}.avg;
%         timelock.var(begchan:endchan,:) = varargin{i}.var;
%         timelock.dof(begchan:endchan,:) = varargin{i}.dof;
%         timelock.label = [timelock.label; varargin{i}.label(:)];
%       end % for varargin
%
%     end
%
%   otherwise
%     error('it is not allowed to concatenate across dimension %s', cfg.appenddim);
% end

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
  opto = cell(1,Ndata);
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
    if hasgrad, timelock.grad = ft_appendsens([], grad{~cellfun(@isempty, grad)}); end
    if haselec, timelock.elec = ft_appendsens([], elec{~cellfun(@isempty, elec)}); end
    if hasopto, timelock.opto = ft_appendsens([], opto{~cellfun(@isempty, opto)}); end
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
    if hasgrad && ~removegrad, timelock.grad = grad{1}; end
    if haselec && ~removeelec, timelock.elec = elec{1}; end
    if hasopto && ~removeopto, timelock.opto = opto{1}; end
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock
