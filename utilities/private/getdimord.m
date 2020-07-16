function dimord = getdimord(data, field, varargin)

% GETDIMORD determine the dimensions and order of a data field in a FieldTrip
% structure.
%
% Use as
%   dimord = getdimord(data, field)
%
% See also GETDIMSIZ, GETDATFIELD, FIXDIMORD

% Copyright (C) 2014-2019, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please note that this function is called from many other FT functions. To avoid
% unwanted recursion, you should avoid (where possible) calling other FT functions
% inside this one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(data, field) && isfield(data, 'avg') && isfield(data.avg, field)
  field = ['avg.' field];
elseif ~isfield(data, field) && isfield(data, 'trial') && isfield(data.trial, field)
  field = ['trial.' field];
elseif ~isfield(data, field)
  ft_error('field "%s" not present in data', field);
end

if strncmp(field, 'avg.', 4)
  prefix = '';
  field = field(5:end); % strip the avg
  data.(field) = data.avg.(field); % copy the avg into the main structure
  data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
  prefix = '(rpt)_';
  field = field(7:end); % strip the trial
  data.(field) = data.trial(1).(field); % copy the first trial into the main structure
  data = rmfield(data, 'trial');
else
  prefix = '';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTEMPT 1: the specific dimord is simply present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, [field 'dimord'])
  dimord = data.([field 'dimord']);
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if not present, we need some additional information about the data strucure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nan means that the value is not known and might remain unknown
% inf means that the value is not known but should be known
ntime     = inf;
nfreq     = inf;
nchan     = inf;
nchancmb  = inf;
nsubj     = nan;
nrpt      = nan;
nrpttap   = nan;
ntopochan = inf;
nspike    = inf; % this is only for the first spike channel
nlag      = nan;
npos      = inf;
nori      = nan; % this will be 3 in many cases, 1 after projectmom, and can be >3 for parcels
ntri      = nan;
ntet      = nan;
nhex      = nan;
ndim1     = nan;
ndim2     = nan;
ndim3     = nan;

% use an anonymous function
assign = @(var, val) assignin('caller', var, val);
% it is possible to pass additional ATTEMPTs such as nrpt, nrpttap, etc
for i=1:2:length(varargin)
  assign(varargin{i}, varargin{i+1});
end

% try to determine the size of each possible dimension in the data
if isfield(data, 'label')
  nchan = length(data.label);
end

if isfield(data, 'labelcmb')
  nchancmb = size(data.labelcmb, 1);
end

if isfield(data, 'time')
  if iscell(data.time) && ~isempty(data.time)
    tmp   = getdimsiz(data, 'time');
    ntime = tmp(3); % raw data may contain variable length trials
  else
    ntime = length(data.time);
  end
end

if isfield(data, 'freq')
  nfreq = length(data.freq);
end

if isfield(data, 'trial') && iscell(data.trial)
  % raw data
  nrpt = length(data.trial);
end

if isfield(data, 'trialtime') && isfield(data, 'timestamp') && isfield(data, 'label')
  % spike data
  nrpt = size(data.trialtime,1);
end

if isfield(data, 'cumtapcnt')
  nrpt = size(data.cumtapcnt,1);
  if numel(data.cumtapcnt)==length(data.cumtapcnt)
    % it is a vector, hence it only represents repetitions
    nrpttap = sum(data.cumtapcnt);
  else
    % it is a matrix, hence it is repetitions by frequencies
    % this happens after  mtmconvol with keeptrials
    nrpttap = sum(data.cumtapcnt,2);
    if any(nrpttap~=nrpttap(1))
      ft_warning('unexpected variation of the number of tapers over trials')
      nrpttap = nan;
    else
      nrpttap = nrpttap(1);
    end
  end
end

if isfield(data, 'pos')
  npos = size(data.pos,1);
elseif isfield(data, 'dim')
  npos = prod(data.dim);
elseif isfield(data, 'leadfield')
  npos = numel(data.leadfield);
elseif isfield(data, 'filter')
  npos = numel(data.filter);
elseif isfield(data, 'inside')
  npos = numel(data.inside);
end

if isfield(data, 'tri')
  ntri = size(data.tri,1);
end

if isfield(data, 'tet')
  ntet = size(data.tet,1);
end

if isfield(data, 'hex')
  nhex = size(data.hex,1);
end

if isfield(data, 'dim')
  ndim1 = data.dim(1);
  ndim2 = data.dim(2);
  ndim3 = data.dim(3);
end

if isfield(data, 'csdlabel')
  % this is used in PCC beamformers
  if length(data.csdlabel)==npos
    % each position has its own labels
    len = cellfun(@numel, data.csdlabel);
    len = len(len~=0);
    if all(len==len(1))
      % they all have the same length
      nori = len(1);
    end
  else
    % one list of labels for all positions
    nori = length(data.csdlabel);
  end
elseif isfield(data, 'mom') && isfield(data, 'inside') && iscell(data.mom)
    % this is used in LCMV beamformers
    size1 = @(x) size(x, 1);
    len = cellfun(size1, data.mom(data.inside));
    if all(len==len(1))
      % they all have the same length
      nori = len(1);
    end
else
  % assume that there are three dipole orientations per source
  nori = 3;
end

if isfield(data, 'topolabel')
  % this is used in ICA and PCA decompositions
  ntopochan = length(data.topolabel);
end

if isfield(data, 'timestamp') && iscell(data.timestamp)
  nspike = length(data.timestamp{1}); % spike data: only for the first channel
end

if isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'lag')) && isfield(data, 'coeffs')
  % mvar data
  nlag = size(data.coeffs,3);
end

% determine the size of the actual data
datsiz = getdimsiz(data, field);

tok = {'subj' 'rpt' 'rpttap' 'chan' 'chancmb' 'freq' 'time' 'topochan' 'lag' 'pos' 'ori' 'tri' 'tet' 'hex' 'dim1' 'dim2' 'dim3'};
siz = [nsubj nrpt nrpttap nchan nchancmb nfreq ntime ntopochan nlag npos nori ntri ntet nhex ndim1 ndim2 ndim3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTEMPT 2: a general dimord is present and might apply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'dimord')
  dimtok  = tokenize(data.dimord, '_');
  if length(dimtok)>length(datsiz) && check_trailingdimsunitlength(data, dimtok((length(datsiz)+1):end))
    % add the trailing singleton dimensions to datsiz, if needed
    datsiz  = [datsiz ones(1,max(0,length(dimtok)-length(datsiz)))];
  end
  if length(dimtok)==length(datsiz) || (length(dimtok)==(length(datsiz)-1) && datsiz(end)==1)
    success = false(size(dimtok));
    for i=1:length(dimtok)
      sel = strcmp(tok, dimtok{i});
      if any(sel) && datsiz(i)==siz(sel)
        success(i) = true;
      elseif strcmp(dimtok{i}, 'subj')
        % the number of subjects cannot be determined, and will be indicated as nan
        success(i) = true;
      elseif strcmp(dimtok{i}, 'rpt')
        % the number of trials is hard to determine, and might be indicated as nan
        success(i) = true;
      end
    end % for
    if all(success)
      dimord = data.dimord;
      return
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTEMPT 3: look at the size of some common fields that are known
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch field
  % the logic for this code is to first check whether the size of a field
  % has an exact match to a potential dimensionality, if not, check for a
  % partial match (ignoring nans)
  
  % note that the case for a cell dimension (typically pos) is handled at
  % the end of this section
  
  case {'pos'}
    if isequalwithoutnans(datsiz, [npos 3])
      dimord = 'pos_unknown';
    end
    
  case {'tri'}
    if isequalwithoutnans(datsiz, [ntri 3])
      dimord = 'tri_unknown';
    end
    
  case {'tet'}
    if isequalwithoutnans(datsiz, [ntet 4])
      dimord = 'tet_unknown';
    end
    
  case {'hex'}
    if isequalwithoutnans(datsiz, [nhex 8])
      dimord = 'hex_unknown';
    end
    
  case {'individual'}
    if isequalwithoutnans(datsiz, [nsubj nchan ntime])
      dimord = 'subj_chan_time';
    end
    
  case {'avg' 'var' 'dof'}
    if isequal(datsiz, [nrpt nchan ntime])
      dimord = 'rpt_chan_time';
    elseif isequal(datsiz, [nchan ntime])
      dimord = 'chan_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan ntime])
      dimord = 'rpt_chan_time';
    elseif isequalwithoutnans(datsiz, [nchan ntime])
      dimord = 'chan_time';
    end
    
  case {'powspctrm' 'fourierspctrm'}
    if isequal(datsiz, [nrpt nchan nfreq ntime])
      dimord = 'rpt_chan_freq_time';
    elseif isequal(datsiz, [nrpt nchan nfreq])
      dimord = 'rpt_chan_freq';
    elseif isequal(datsiz, [nchan nfreq ntime])
      dimord = 'chan_freq_time';
    elseif isequal(datsiz, [nchan nfreq])
      dimord = 'chan_freq';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nfreq ntime])
      dimord = 'rpt_chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nfreq])
      dimord = 'rpt_chan_freq';
    elseif isequalwithoutnans(datsiz, [nchan nfreq ntime])
      dimord = 'chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nchan nfreq])
      dimord = 'chan_freq';
    end
    
  case {'crsspctrm' 'cohspctrm'}
    if isequal(datsiz, [nrpt nchancmb nfreq ntime])
      dimord = 'rpt_chancmb_freq_time';
    elseif isequal(datsiz, [nrpt nchancmb nfreq])
      dimord = 'rpt_chancmb_freq';
    elseif isequal(datsiz, [nchancmb nfreq ntime])
      dimord = 'chancmb_freq_time';
    elseif isequal(datsiz, [nchancmb nfreq])
      dimord = 'chancmb_freq';
    elseif isequal(datsiz, [nrpt nchan nchan nfreq ntime])
      dimord = 'rpt_chan_chan_freq_time';
    elseif isequal(datsiz, [nrpt nchan nchan nfreq])
      dimord = 'rpt_chan_chan_freq';
    elseif isequal(datsiz, [nchan nchan nfreq ntime])
      dimord = 'chan_chan_freq_time';
    elseif isequal(datsiz, [nchan nchan nfreq])
      dimord = 'chan_chan_freq';
    elseif isequal(datsiz, [npos nori])
      dimord = 'pos_ori';
    elseif isequal(datsiz, [npos 1])
      dimord = 'pos';
    elseif isequalwithoutnans(datsiz, [nrpt nchancmb nfreq ntime])
      dimord = 'rpt_chancmb_freq_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchancmb nfreq])
      dimord = 'rpt_chancmb_freq';
    elseif isequalwithoutnans(datsiz, [nchancmb nfreq ntime])
      dimord = 'chancmb_freq_time';
    elseif isequalwithoutnans(datsiz, [nchancmb nfreq])
      dimord = 'chancmb_freq';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nchan nfreq ntime])
      dimord = 'rpt_chan_chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nchan nfreq])
      dimord = 'rpt_chan_chan_freq';
    elseif isequalwithoutnans(datsiz, [nchan nchan nfreq ntime])
      dimord = 'chan_chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nchan nchan nfreq])
      dimord = 'chan_chan_freq';
    elseif isequalwithoutnans(datsiz, [npos nori])
      dimord = 'pos_ori';
    elseif isequalwithoutnans(datsiz, [npos 1])
      dimord = 'pos';
    end
    
  case {'cov' 'coh' 'csd' 'noisecov' 'noisecsd'}
    % these occur in timelock and in source structures
    if isequal(datsiz, [nrpt nchan nchan])
      dimord = 'rpt_chan_chan';
    elseif isequal(datsiz, [nchan nchan])
      dimord = 'chan_chan';
    elseif isequal(datsiz, [npos nori nori])
      dimord = 'pos_ori_ori';
    elseif isequal(datsiz, [npos nrpt nori nori])
      dimord = 'pos_rpt_ori_ori';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nchan])
      dimord = 'rpt_chan_chan';
    elseif isequalwithoutnans(datsiz, [nchan nchan])
      dimord = 'chan_chan';
    elseif isequalwithoutnans(datsiz, [npos nori nori])
      dimord = 'pos_ori_ori';
    elseif isequalwithoutnans(datsiz, [npos nrpt nori nori])
      dimord = 'pos_rpt_ori_ori';
    end
    
  case {'tf'}
    if isequal(datsiz, [npos nfreq ntime])
      dimord = 'pos_freq_time';
    end
    
  case {'pow' 'noise' 'rv' 'nai' 'kurtosis'}
    if isequal(datsiz, [npos ntime])
      dimord = 'pos_time';
    elseif isequal(datsiz, [npos nfreq])
      dimord = 'pos_freq';
    elseif isequal(datsiz, [npos nrpt])
      dimord = 'pos_rpt';
    elseif isequal(datsiz, [nrpt npos ntime])
      dimord = 'rpt_pos_time';
    elseif isequal(datsiz, [nrpt npos nfreq])
      dimord = 'rpt_pos_freq';
    elseif isequal(datsiz, [npos 1]) % in case there are no repetitions
      dimord = 'pos';
    elseif isequalwithoutnans(datsiz, [npos ntime])
      dimord = 'pos_time';
    elseif isequalwithoutnans(datsiz, [npos nfreq])
      dimord = 'pos_freq';
    elseif isequalwithoutnans(datsiz, [npos nrpt])
      dimord = 'pos_rpt';
    elseif isequalwithoutnans(datsiz, [nrpt npos ntime])
      dimord = 'rpt_pos_time';
    elseif isequalwithoutnans(datsiz, [nrpt npos nfreq])
      dimord = 'rpt_pos_freq';
    end
    
  case {'mom' 'itc' 'aa' 'stat','pval' 'statitc' 'pitc'}
    if isequal(datsiz, [npos nori nrpt])
      dimord = 'pos_ori_rpt';
    elseif isequal(datsiz, [npos nori ntime])
      dimord = 'pos_ori_time';
    elseif isequal(datsiz, [npos nori nfreq])
      dimord = 'pos_ori_nfreq';
    elseif isequal(datsiz, [npos ntime])
      dimord = 'pos_time';
    elseif isequal(datsiz, [npos nfreq])
      dimord = 'pos_freq';
    elseif isequal(datsiz, [npos 3])
      dimord = 'pos_ori';
    elseif isequal(datsiz, [npos 1])
      dimord = 'pos';
    elseif isequal(datsiz, [npos nrpt])
      dimord = 'pos_rpt';
    elseif isequalwithoutnans(datsiz, [npos nori nrpt])
      dimord = 'pos_ori_rpt';
    elseif isequalwithoutnans(datsiz, [npos nori nrpttap])
      dimord = 'pos_ori_rpttap';
    elseif isequalwithoutnans(datsiz, [npos nori ntime])
      dimord = 'pos_ori_time';
    elseif isequalwithoutnans(datsiz, [npos nori nfreq])
      dimord = 'pos_ori_nfreq';
    elseif isequalwithoutnans(datsiz, [npos ntime])
      dimord = 'pos_time';
    elseif isequalwithoutnans(datsiz, [npos nfreq])
      dimord = 'pos_freq';
    elseif isequalwithoutnans(datsiz, [npos 3])
      dimord = 'pos_ori';
    elseif isequalwithoutnans(datsiz, [npos 1])
      dimord = 'pos';
    elseif isequalwithoutnans(datsiz, [npos nrpt])
      dimord = 'pos_rpt';
    elseif isequalwithoutnans(datsiz, [npos nrpt nori ntime])
      dimord = 'pos_rpt_ori_time';
    elseif isequalwithoutnans(datsiz, [npos nrpt 1 ntime])
      dimord = 'pos_rpt_ori_time';
    elseif isequal(datsiz, [npos nfreq ntime])
      dimord = 'pos_freq_time';
    end
    
  case {'filter'}
    if isequalwithoutnans(datsiz, [npos nori nchan]) || (isequal(datsiz([1 2]), [npos nori]) && isinf(nchan))
      dimord = 'pos_ori_chan';
    end
    
  case {'leadfield'}
    if isequalwithoutnans(datsiz, [npos nchan nori]) || (isequal(datsiz([1 3]), [npos nori]) && isinf(nchan))
      dimord = 'pos_chan_ori';
    end
    
  case {'ori' 'eta'}
    if isequal(datsiz, [npos nori]) || isequal(datsiz, [npos nori 1]) || isequal(datsiz, [npos 3]) || isequal(datsiz, [npos 3 1])
      dimord = 'pos_ori';
    elseif isequal(datsiz, [npos 1 nori]) || isequal(datsiz, [npos 1 3])
      dimord = 'pos_unknown_ori';
    end
    
  case {'csdlabel'}
    if isequal(datsiz, [npos nori]) || isequal(datsiz, [npos 3])
      dimord = 'pos_ori';
    end
    
  case {'trial'}
    if ~iscell(data.(field)) && isequalwithoutnans(datsiz, [nrpt nchan ntime])
      dimord = 'rpt_chan_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan ntime])
      dimord = '{rpt}_chan_time';
    elseif isequalwithoutnans(datsiz, [nchan nspike]) || isequalwithoutnans(datsiz, [nchan 1 nspike])
      dimord = '{chan}_spike';
    end
    
  case {'sampleinfo' 'trialinfo' 'trialtime'}
    if isequalwithoutnans(datsiz, [nrpt nan])
      dimord = 'rpt_other';
    end
    
  case {'cumtapcnt' 'cumsumcnt'}
    if isequalwithoutnans(datsiz, [nrpt 1])
      dimord = 'rpt';
    elseif isequalwithoutnans(datsiz, [nrpt nfreq])
      dimord = 'rpt_freq';
    elseif isequalwithoutnans(datsiz, [nrpt nan])
      dimord = 'rpt_other';
    end
    
  case {'topo'}
    if isequalwithoutnans(datsiz, [ntopochan nchan])
      dimord = 'topochan_chan';
    end
    
  case {'unmixing'}
    if isequalwithoutnans(datsiz, [nchan ntopochan])
      dimord = 'chan_topochan';
    end
    
  case {'anatomy' 'inside'}
    if isfield(data, 'dim') && isequal(datsiz, data.dim)
      dimord = 'dim1_dim2_dim3';
    elseif isequalwithoutnans(datsiz, [npos 1]) || isequalwithoutnans(datsiz, [1 npos])
      dimord = 'pos';
    end
    
  case {'timestamp'}
    if iscell(data.(field)) && isfield(data, 'label') && datsiz(1)==nchan
      dimord = '{chan}_spike';
    end
    
  case {'time'}
    if iscell(data.(field)) && isfield(data, 'label') && datsiz(1)==nrpt
      dimord = '{rpt}_time';
    elseif isvector(data.(field)) && isequal(datsiz, [1 ntime ones(1,numel(datsiz)-2)])
      dimord = 'time';
    elseif iscell(data.(field)) && isfield(data, 'label') && isfield(data, 'timestamp') && isequal(getdimsiz(data, 'timestamp'), datsiz) && datsiz(1)==nchan
      dimord = '{chan}_spike';
    end
    
  case {'freq'}
    if iscell(data.(field)) && isfield(data, 'label') && datsiz(1)==nrpt
      dimord = '{rpt}_freq';
    elseif isvector(data.(field)) && isequal(datsiz, [1 nfreq ones(1,numel(datsiz)-2)])
      dimord = 'freq';
    end
    
  case {'chantype', 'chanunit'}
    if numel(data.(field))==nchan
      dimord = 'chan';
    end
    
  otherwise
    if isfield(data, 'dim') && isequal(datsiz, data.dim)
      dimord = 'dim1_dim2_dim3';
    end
    
end % switch field

% deal with possible first pos which is a cell
if exist('dimord', 'var') && strcmp(dimord(1:3), 'pos') && iscell(data.(field))
  dimord = ['{pos}' dimord(4:end)];
end

if ~exist('dimord', 'var')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ATTEMPT 4: there is only one way that the dimensions can be interpreted
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dimtok = cell(size(datsiz));
  
  for i=1:length(datsiz)
    sel = find(siz==datsiz(i));
    if length(sel)==1
      % there is exactly one corresponding dimension
      dimtok{i} = tok{sel};
    else
      % there are zero or multiple corresponding dimensions
      dimtok{i} = [];
    end
  end
  
  if all(~cellfun(@isempty, dimtok))
    % each of the dimensions matches uniquely with a single known size
    if iscell(data.(field))
      dimtok{1} = ['{' dimtok{1} '}'];
    end
    dimord = sprintf('%s_', dimtok{:});
    dimord = dimord(1:end-1);
    return
  elseif ~isempty(dimtok{1}) && numel(datsiz)==2 && datsiz(2)==1
    % it is often impossible to determine or specify what the 2nd dimension is in a Nx1 matrix
    dimord = dimtok{1};
    return
  end
end % if dimord does not exist

if ~exist('dimord', 'var')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ATTEMPT 5: compare the size with the known size of each dimension
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sel = ~isnan(siz) & ~isinf(siz);
  % nan means that the value is not known and might remain unknown
  % inf means that the value is not known and but should be known
  if length(unique(siz(sel)))==length(siz(sel))
    % this should only be done if there is no chance of confusing dimensions
    dimtok = cell(size(datsiz));
    dimtok(datsiz==npos)      = {'pos'};
    dimtok(datsiz==nori)      = {'ori'};
    dimtok(datsiz==nrpttap)   = {'rpttap'};
    dimtok(datsiz==nrpt)      = {'rpt'};
    dimtok(datsiz==nsubj)     = {'subj'};
    dimtok(datsiz==nchancmb)  = {'chancmb'};
    dimtok(datsiz==nchan)     = {'chan'};
    dimtok(datsiz==nfreq)     = {'freq'};
    dimtok(datsiz==ntime)     = {'time'};
    dimtok(datsiz==ndim1)     = {'dim1'};
    dimtok(datsiz==ndim2)     = {'dim2'};
    dimtok(datsiz==ndim3)     = {'dim3'};
    
    if isempty(dimtok{end}) && datsiz(end)==1
      % remove the unknown trailing singleton dimension
      dimtok = dimtok(1:end-1);
    elseif isequal(dimtok{1}, 'pos') && isempty(dimtok{2}) && datsiz(2)==1
      % remove the unknown leading singleton dimension
      dimtok(2) = [];
    end
    
    if all(~cellfun(@isempty, dimtok))
      if iscell(data.(field))
        dimtok{1} = ['{' dimtok{1} '}'];
      end
      dimord = sprintf('%s_', dimtok{:});
      dimord = dimord(1:end-1);
      return
    end
  end
end % if dimord does not exist

if ~exist('dimord', 'var')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ATTEMPT 6: check whether it is a 3-D volume
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isequal(datsiz, [ndim1 ndim2 ndim3])
    dimord = 'dim1_dim2_dim3';
    return
  elseif isfield(data, 'pos') && prod(datsiz)==size(data.pos, 1)
    dimord = 'dim1_dim2_dim3';
    return
  end
end % if dimord does not exist



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL RESORT: return "unknown" for all unknown dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('dimord', 'var')
  % this should not happen
  % if it does, it might help in diagnosis to have a very informative warning message
  % since there have been problems with trials not being selected correctly due to the warning going unnoticed
  % it is better to throw an error than a warning
  warning_dimord_could_not_be_determined(field, data);
  
  dimtok(cellfun(@isempty, dimtok)) = {'unknown'};
  if all(~cellfun(@isempty, dimtok))
    if iscell(data.(field))
      dimtok{1} = ['{' dimtok{1} '}'];
    end
    dimord = sprintf('%s_', dimtok{:});
    dimord = dimord(1:end-1);
  end
end

% add '(rpt)' in case of source.trial
dimord = [prefix dimord];


end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warning_dimord_could_not_be_determined(field,data)
  msg=sprintf('could not determine dimord of "%s" in:',field);

  if isempty(which('evalc'))
    % May not be available in Octave
    content=sprintf('object of type ''%s''',class(data));
  else
    % in Octave, disp typically shows full data arrays which can result in
    % very long output. Here we take out the middle part of the output if
    % the output is very long (more than 40 lines)
    full_content=evalc('disp(data)');
    max_pre_post_lines=20;

    newline_pos=find(full_content==newline);
    newline_pos=newline_pos(max_pre_post_lines:(end-max_pre_post_lines));

    if numel(newline_pos)>=2
      pre_end=newline_pos(1)-1;
      post_end=newline_pos(end)+1;

      content=sprintf('%s\n\n... long output omitted ...\n\n%s',...
                                full_content(1:pre_end),...
                                full_content(post_end:end));
    else
      content=full_content;
    end
  end

  msg = sprintf('%s\n\n%s', msg, content);
  ft_warning(msg);
end % function warning_dimord_could_not_be_determined


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = isequalwithoutnans(a, b)

% this is *only* used to compare matrix sizes, so we can ignore any singleton last dimension
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
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_trailingdimsunitlength(data, dimtok)

ok = false;
for k = 1:numel(dimtok)
  switch dimtok{k}
    case 'chan'
      ok = numel(data.label)==1;
    otherwise
      if isfield(data, dimtok{k}) % check whether field exists
        ok = numel(data.(dimtok{k}))==1;
      end
  end
  if ok
    break;
  end
end

end % function check_trailingdimsunitlength
