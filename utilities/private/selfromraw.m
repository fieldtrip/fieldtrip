function [data] = selfromraw(data, varargin)

% FIXME this function is not documented

selrpt  = ft_getopt(varargin, 'rpt');
selchan = ft_getopt(varargin, 'chan');
seltim  = ft_getopt(varargin, 'latency');

if ~isempty(selrpt)
  fprintf('selecting %d trials\n', numel(selrpt));
  data.trial  = data.trial(selrpt);
  data.time   = data.time(selrpt);
  
  % this is new style keeping track of sampleinfo and trialinfo
  if isfield(data, 'sampleinfo'), data.sampleinfo = data.sampleinfo(selrpt, :); end
  if isfield(data, 'trialinfo'),  data.trialinfo  = data.trialinfo(selrpt,  :); end
  if isfield(data, 'offset'),     data.offset     = data.offset(selrpt);        end

end

if ~isempty(selchan)
  for k = 1:numel(data.trial)
    data.trial{k} = data.trial{k}(selchan, :);
  end
  data.label = data.label(selchan);
end

if ~isempty(seltim)
  
  oktrial = true(numel(data.trial), 1);
  for k = 1:numel(data.trial)
    ok = data.time{k}>=seltim(1) & data.time{k}<=seltim(2);
    
    if sum(ok)>0
      % there are samples in trial k that fall within latency
      
      data.trial{k} = data.trial{k}(:,ok);
      data.time{k}  = data.time{k}(ok);
      
      % how many samples at the beginning and/or end were removed?
      beginLatency = find(ok,1,'first') - 1;
      endLatency = numel(ok) - find(ok,1,'last');
      
      % update sampleinfo to reflect this
      data.sampleinfo(k,1) = data.sampleinfo(k,1) + beginLatency;
      data.sampleinfo(k,2) = data.sampleinfo(k,2) - endLatency;
      
    else
      
      % whole trial falls outside latency
      oktrial(k) = false;
      
    end
  end
  
  % remove all trials that had no samples within latency
  data = selfromraw(data, 'rpt', find(oktrial));
  
end
