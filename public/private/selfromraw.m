function [data] = selfromraw(data, varargin)

selrpt  = keyval('rpt',  varargin{:}); selectrpt  = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'rpt'));
selchan = keyval('chan', varargin{:}); selectchan = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'chan'));

if selectrpt,
  fprintf('selecting %d trials\n', numel(selrpt));
  data.trial  = data.trial(selrpt);
  data.time   = data.time(selrpt);
  
  % this is new style keeping track of trialdef and trialinfo
  if isfield(data, 'trialdef'),  data.trialdef  = data.trialdef(selrpt,  :); end
  if isfield(data, 'trialinfo'), data.trialinfo = data.trialinfo(selrpt, :); end
  if isfield(data, 'offset'),    data.offset    = data.offset(selrpt);       end

  % this is old style keeping track of trl-matrix
  if isfield(data, 'cfg')
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end

  if isempty(trl)
    warning('could not locate the correct trial definition ''trl'' in the data structure');
  else
    data.cfg.trlold = trl;
    data.cfg.trl    = trl(selrpt,:); 
  end
  
end

if selectchan,
  for k = 1:numel(data.trial)
    data.trial{k} = data.trial{k}(selchan, :);
  end
  data.label = data.label(selchan);
end
