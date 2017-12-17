function data = append_common(cfg, varargin)

% APPEND_COMMON is used for concatenating raw, timelock or freq data
%
% The general bookkeeping and the correct specification of the cfg
% should be taken care of by the calling function.
%
% See FT_APPENDDATA, T_APPENDTIMELOCK, FT_APPENDFREQ

% Copyright (C) 2017, Robert Oostenveld
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

% when present, these must be consistent in all inputs
hastime       = isfield(varargin{1}, 'time');
hasfreq       = isfield(varargin{1}, 'freq');
hastrialinfo  = isfield(varargin{1}, 'trialinfo');
hassampleinfo = isfield(varargin{1}, 'sampleinfo');
hascumsumcnt  = isfield(varargin{1}, 'cumsumcnt');
hascumtapcnt  = isfield(varargin{1}, 'cumtapcnt');
hastopolabel  = isfield(varargin{1}, 'topolabel');
hastopo       = isfield(varargin{1}, 'topo');
hasunmixing   = isfield(varargin{1}, 'unmixing');
for i=2:numel(varargin)
  hastime       = hastime       && isfield(varargin{i}, 'time');
  hasfreq       = hasfreq       && isfield(varargin{i}, 'freq');
  hastrialinfo  = hastrialinfo  && isfield(varargin{i}, 'trialinfo');
  hassampleinfo = hassampleinfo && isfield(varargin{i}, 'sampleinfo');
  hascumsumcnt  = hascumsumcnt  && isfield(varargin{i}, 'cumsumcnt');
  hascumtapcnt  = hascumtapcnt  && isfield(varargin{i}, 'cumtapcnt');
  hastopolabel  = hastopolabel  && isfield(varargin{i}, 'topolabel');
  hastopo       = hastopo       && isfield(varargin{i}, 'topo');
  hasunmixing   = hasunmixing   && isfield(varargin{i}, 'unmixing');
end

% these can be present in a subset of the inputs, e.g. for appending EEG and MEG data
hasgrad = false;
haselec = false;
hasopto = false;
for i=1:length(varargin)
  hasgrad = hasgrad || isfield(varargin{i}, 'grad');
  haselec = haselec || isfield(varargin{i}, 'elec');
  hasopto = hasopto || isfield(varargin{i}, 'opto');
end

if hastopolabel || hastopo || hasunmixing
  identical = true;
  for i=2:numel(varargin)
    if hastopolabel, identical = identical && isequal(varargin{1}.topolabel, varargin{i}.topolabel); end
    if hastopo,      identical = identical && isequal(varargin{1}.topo,      varargin{i}.topo);      end
    if hasunmixing,  identical = identical && isequal(varargin{1}.unmixing,  varargin{i}.unmixing);  end
  end
  if strcmp(cfg.appenddim, 'chan')
    % It is possible to combine the component timeseries of different decompositions
    % in the same dataset. In principle this could be improved by also concatenating
    % the topo and unmixing along the correct dimension. However, at the moment the
    % topo/unmixing are discarded.
    ft_warning('discarding ICA/PCA topographies and/or unmixing matrix');
  else
    % only proceed if the ICA/PCA topographies and/or unmixing matrix is identical in all datasets
    assert(identical, 'cannot append data from different ICA/PCA decompositions');
  end
end

switch cfg.appenddim
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'chan'
    assert(checkchan(varargin{:}, 'unique'), 'not all channels are unique');
    % remember the original channel labels in each input
    oldlabel = cell(size(varargin));
    for i=1:numel(varargin)
      oldlabel{i} =  varargin{i}.label;
    end
    
    % determine the union of all input data
    tmpcfg = keepfields(cfg, {'tolerance', 'channel', 'showcallinfo'});
    tmpcfg.select = 'union';
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the union of all input data
    data = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord'});
    
    % keep these fields (when identical)
    fn = {'trialinfo' 'sampleinfo', 'cumsumcnt', 'cumtapcnt'};
    for i=1:numel(fn)
      keepfield = isfield(varargin{1}, fn{i});
      for j=1:numel(varargin)
        if ~isfield(varargin{j}, fn{i}) || ~isequaln(varargin{j}.(fn{i}), varargin{1}.(fn{i}))
          keepfield = false;
          break
        end
      end
      if keepfield
        data.(fn{i}) = varargin{1}.(fn{i});
      end
    end % for each of the fields to keep
    
    for i=1:numel(cfg.parameter)
      dimsiz = getdimsiz(varargin{1}, cfg.parameter{i});
      switch getdimord(varargin{1}, cfg.parameter{i})
        case {'chan_chan'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            data.(cfg.parameter{i})(chansel,chansel) = varargin{j}.(cfg.parameter{i})(chansel,chansel);
          end
          
        case {'rpt_chan_chan'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            data.(cfg.parameter{i})(:,chansel,chansel) = varargin{j}.(cfg.parameter{i})(:,chansel,chansel);
          end
          
        case {'chan' 'chan_time' 'chan_freq'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            data.(cfg.parameter{i})(chansel,:) = varargin{j}.(cfg.parameter{i})(chansel,:);
          end
          
        case {'rpt_chan_time' 'subj_chan_time' 'rpt_chan_freq' 'subj_chan_freq'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            chansel = match_str(varargin{j}.label, oldlabel{j});
            data.(cfg.parameter{i})(:,chansel,:) = varargin{j}.(cfg.parameter{i})(:,chansel,:);
          end
          
        otherwise
          % do not concatenate this field
      end % switch
    end % for cfg.parameter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case {'time' 'freq'}
    
    % remember the original axes in each input
    if hasfreq
      oldfreq = cell(size(varargin));
      for i=1:numel(varargin)
        oldfreq{i} = varargin{i}.freq;
      end
    end
    if hastime
      oldtime = cell(size(varargin));
      for i=1:numel(varargin)
        oldtime{i} =  varargin{i}.time;
      end
    end
    
    % determine the union of all input data
    tmpcfg = keepfields(cfg, {'tolerance', 'channel'});
    tmpcfg.select = 'union';
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the union of all input data
    data = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord', 'topo', 'unmixing', 'topolabel'});
    
    % keep the trialinfo (when identical)
    % note that we are NOT keeing the sampleinfo, cumsumcnt, cumtapcnt
    fn = {'trialinfo'};
    for i=1:numel(fn)
      keepfield = isfield(varargin{1}, fn{i});
      for j=1:numel(varargin)
        if ~isfield(varargin{j}, fn{i}) || ~isequal(varargin{j}.(fn{i}), varargin{1}.(fn{i}))
          keepfield = false;
          break
        end
      end
      if keepfield
        data.(fn{i}) = varargin{1}.(fn{i});
      end
    end % for each of the fields to keep
    
    for i=1:numel(cfg.parameter)
      dimsiz = getdimsiz(varargin{1}, cfg.parameter{i});
      switch getdimord(varargin{1}, cfg.parameter{i})
        case 'chan_time'
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timesel = match_val(varargin{j}.time, oldtime{j});
            data.(cfg.parameter{i})(:,timesel) = varargin{j}.(cfg.parameter{i})(:,timesel);
          end
          
        case 'chan_freq'
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            data.(cfg.parameter{i})(:,freqsel) = varargin{j}.(cfg.parameter{i})(:,freqsel);
          end
          
        case 'chan_freq_time'
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            timesel = match_val(varargin{j}.time, oldtime{j});
            data.(cfg.parameter{i})(:,freqsel,timesel) = varargin{j}.(cfg.parameter{i})(:,freqsel,timesel);
          end
          
        case {'rpt_chan_time' 'subj_chan_time'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            timesel = match_val(varargin{j}.time, oldtime{j});
            data.(cfg.parameter{i})(:,:,timesel) = varargin{j}.(cfg.parameter{i})(:,:,timesel);
          end
          
        case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            data.(cfg.parameter{i})(:,:,freqsel) = varargin{j}.(cfg.parameter{i})(:,:,freqsel);
          end
          
        case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
          data.(cfg.parameter{i}) = nan(dimsiz);
          for j=1:numel(varargin)
            freqsel = match_val(varargin{j}.freq, oldfreq{j});
            timesel = match_val(varargin{j}.time, oldtime{j});
            data.(cfg.parameter{i})(:,:,freqsel,timesel) = varargin{j}.(cfg.parameter{i})(:,:,freqsel,timesel);
          end
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for cfg.parameter
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'rpt'
    
    % determine the intersection of all input data
    tmpcfg = keepfields(cfg, {'tolerance', 'channel'});
    tmpcfg.select = 'intersect';
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:numel(varargin)
      [cfg, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    
    % start with the intersection of all input data
    data = keepfields(varargin{1}, {'label', 'time', 'freq', 'dimord', 'topo', 'unmixing', 'topolabel'});
    if numel(cfg.parameter)>0
      % this check should not be done if there is no data to append, this happens when called from ft_appenddata
      assert(numel(data.label)>0);
    end
    if hastime, assert(numel(data.time)>0); end
    if hasfreq, assert(numel(data.freq)>0); end
    
    % also append these when present
    if hastrialinfo,  cfg.parameter{end+1} = 'trialinfo';  end
    if hassampleinfo, cfg.parameter{end+1} = 'sampleinfo'; end
    if hascumsumcnt,  cfg.parameter{end+1} = 'cumsumcnt';  end
    if hascumtapcnt,  cfg.parameter{end+1} = 'cumtapcnt';  end
    
    for i=1:numel(cfg.parameter)
      dimsiz = getdimsiz(varargin{1}, cfg.parameter{i});
      switch getdimord(varargin{1}, cfg.parameter{i})
        case {'chan' 'chan_time' 'chan_freq' 'chan_chan' 'chan_freq_time' 'chan_chan_freq' 'chan_chan_time' 'chan_chan_freq_time'}
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            % add a singleton dimension to the beginning
            dat{j} = reshape(varargin{j}.(cfg.parameter{i}), [1, dimsiz]);
          end
          data.(cfg.parameter{i}) = cat(1, dat{:});
          
        case {'rpt' 'rpt_chan' 'rpt_chan_time' 'rpt_chan_freq' 'rpt_chan_chan' 'rpt_chan_freq_time' 'rpttap_chan_freq' 'rpttap_chan_freq_time' 'rpt_other'}
          dat = cell(size(varargin));
          for j=1:numel(varargin)
            dat{j} = varargin{j}.(cfg.parameter{i});
          end
          data.(cfg.parameter{i}) = cat(1, dat{:});
          
        otherwise
          % do not concatenate this field
          
      end % switch
    end % for cfg.parameter
    
  otherwise
    ft_error('unsupported cfg.appenddim');
end

if isfield(data, 'dimord')
  dimtok = tokenize(data.dimord, '_');
  if strcmp(cfg.appenddim, 'rpt') && ~any(strcmp(dimtok{1}, {'rpt', 'rpttap', 'subj'}))
    data.dimord = ['rpt_' data.dimord];
  end
end

if hasgrad || haselec || hasopto
  % gather the sensor definitions from all inputs
  grad = cell(size(varargin));
  elec = cell(size(varargin));
  opto = cell(size(varargin));
  for j=1:length(varargin)
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
    if hasgrad, data.grad = ft_appendsens([], grad{~cellfun(@isempty, grad)}); end
    if haselec, data.elec = ft_appendsens([], elec{~cellfun(@isempty, elec)}); end
    if hasopto, data.opto = ft_appendsens([], opto{~cellfun(@isempty, opto)}); end
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
    if hasgrad && ~removegrad, data.grad = grad{1}; end
    if haselec && ~removeelec, data.elec = elec{1}; end
    if hasopto && ~removeopto, data.opto = opto{1}; end
  end
end
