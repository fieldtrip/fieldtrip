function dimord = getdimord(data, field)

% GETDIMORD
%
% Use as
%   dimord = getdimord(data, field)
%
% See also GETDIMSIZ

if ~isfield(data, field)
  error('field "%s" not present in data', field);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 1: the specific dimord is simply present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, [field 'dimord'])
  dimord = data.([field 'dimord']);
  return
end

% if not present, we need some additional information about the data strucure

% nan means that the value is not known and might remain unknown
% inf means that the value is not known and but should be known
ntime     = inf;
nfreq     = inf;
nchan     = inf;
nchancmb  = inf;
nrpt      = nan; 
nrpttap   = nan;
npos      = inf;
ntopochan = inf;

% try to determine the size of each possible dimension in the data
if isfield(data, 'label')
  nchan = length(data.label);
end
if isfield(data, 'labelcmb')
  nchancmb = size(data.labelcmb, 1);
end
if isfield(data, 'time')
  ntime = length(data.time);
end
if isfield(data, 'freq')
  nfreq = length(data.freq);
end
if isfield(data, 'trial') && iscell(data.trial)
  nrpt = length(data.trial);
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
      warning('unexpected variation of the number of tapers over trials')
      nrpttap = nan;
    else
      nrpttap = nrpttap(1);
    end
  end
end
if isfield(data, 'pos')
  npos = size(data.pos,1);
end
if isfield(data, 'topolabel')
  ntopochan = length(data.topolabel);
end

% determine the size of the actual data
datsiz = getdimsiz(data, field);

tok = {'rpt' 'rpttap' 'chan' 'chancmb' 'freq' 'time' 'pos' 'topochan'};
siz = [nrpt nrpttap nchan nchancmb nfreq ntime npos ntopochan];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 2: a general dimord is present and might apply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'dimord')
  dimtok = cell(size(datsiz));
  
  for i=1:length(datsiz)
    sel = find(siz==datsiz(i));
    if length(sel)==1
      dimtok{i} = tok{sel};
    else
      dimtok{i} = [];
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 3: look at the size of some common fields that are known
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch field
  case {'avg' 'var' 'dof'}
    if isequalwithoutnans(datsiz, [nrpt nchan ntime])
      dimord = 'rpt_chan_time';
    elseif isequalwithoutnans(datsiz, [nchan ntime])
      dimord = 'chan_time';
    end
    
  case {'powspctrm' 'fourierspctrm'}
    if isequalwithoutnans(datsiz, [nrpt nchan nfreq ntime])
      dimord = 'rpt_chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan nfreq])
      dimord = 'rpt_chan_freq';
    elseif isequalwithoutnans(datsiz, [nchan nfreq ntime])
      dimord = 'chan_freq_time';
    elseif isequalwithoutnans(datsiz, [nchan nfreq])
      dimord = 'chan_freq';
    end
    
  case {'crsspctrm' 'cohspctrm'}
    if isequalwithoutnans(datsiz, [nrpt nchancmb nfreq ntime])
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
    end
    
  case {'cov'}
    if isequalwithoutnans(datsiz, [nrpt nchan nchan])
      dimord = 'rpt_chan_chan';
    elseif isequalwithoutnans(datsiz, [nchan nchan])
      dimord = 'chan_chan';
    end
    
  case {'trial'}
    if iscell(data.trial)
      dimord = '{rpt}_chan_time';
    elseif isequalwithoutnans(datsiz, [nrpt nchan ntime])
      dimord = 'rpt_chan_time';
    end
    
  case {'sampleinfo' 'trialinfo'}
    if isequalwithoutnans(datsiz, [nrpt nan])
      dimord = 'rpt_unknown';
    end
    
  case {'cumtapcnt' 'cumsumcnt'}
    if isequalwithoutnans(datsiz, [nrpt nan])
      dimord = 'rpt_unknown';
    end
        
  case {'topo'}
    if isequalwithoutnans(datsiz, [ntopochan nchan])
      dimord = 'topochan_chan';
    end
    
  case {'unmixing'}
    if isequalwithoutnans(datsiz, [nchan ntopochan])
      dimord = 'chan_topochan';
    end

  case {'inside'}
    if isequalwithoutnans(datsiz, [npos])
      dimord = 'pos';
    end
    
  otherwise
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTION 4: compare the size with the known size of each dimension
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if length(unique(siz(~isnan(siz))))==length(siz(~isnan(siz)))
      % this should only be done if there is no chance of confusing dimensions
      dimtok = cell(size(datsiz));
      dimtok(datsiz==npos)      = {'pos'};
      dimtok(datsiz==nrpttap)   = {'rpttap'};
      dimtok(datsiz==nrpt)      = {'rpt'};
      dimtok(datsiz==nchancmb)  = {'chancmb'};
      dimtok(datsiz==nchan)     = {'chan'};
      dimtok(datsiz==nfreq)     = {'freq'};
      dimtok(datsiz==ntime)     = {'time'};
      if all(~cellfun(@isempty, dimtok))
        if iscell(data.(field))
          dimtok{1} = ['{' dimtok{1} '}'];
        end
        dimord = sprintf('%s_', dimtok{:});
        dimord = dimord(1:end-1);
        return
      end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 5: return "unknown_unknown"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('dimord', 'var')
  % this should not happen
  % if it does, it might help in diagnosis to have a very informative warning message
  warning('could not determine dimord of "%s" in the following data', field)
  disp(data);
  
  dimtok = repmat({'unknown'}, size(datsiz));
  if all(~cellfun(@isempty, dimtok))
    if iscell(data.(field))
      dimtok{1} = ['{' dimtok{1} '}'];
    end
    dimord = sprintf('%s_', dimtok{:});
    dimord = dimord(1:end-1);
  end
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
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

