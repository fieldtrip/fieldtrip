function [data] = selfromraw(data, varargin)

selrpt  = keyval('rpt',  varargin); selectrpt  = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'rpt'));
selchan = keyval('chan', varargin); selectchan = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'chan'));
seltim  = keyval('latency', varargin);

if selectrpt,
  fprintf('selecting %d trials\n', numel(selrpt));
  data.trial  = data.trial(selrpt);
  data.time   = data.time(selrpt);
  
  % this is new style keeping track of sampleinfo and trialinfo
  if isfield(data, 'sampleinfo'), data.sampleinfo = data.sampleinfo(selrpt, :); end
  if isfield(data, 'trialinfo'),  data.trialinfo  = data.trialinfo(selrpt,  :); end
  if isfield(data, 'offset'),     data.offset     = data.offset(selrpt);        end

end

if selectchan,
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
      data.trial{k} = data.trial{k}(:,ok);
      data.time{k}  = data.time{k}(ok);
    else
      oktrial(k) = false;
    end
  end
  
  % remove sampleinfo if present
  % FIXME change this into update sampleinfo
  if isfield(data, 'sampleinfo'), data = rmfield(data, 'sampleinfo'); end
  data = selfromraw(data, 'rpt', find(oktrial));
end
