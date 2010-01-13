function [data] = selfromraw(data, varargin)

selrpt  = keyval('rpt',  varargin{:});
selchan = keyval('chan', varargin{:});

if ~isempty(selrpt)
  fprintf('selecting %d trials\n', numel(selrpt));
  data.trial  = data.trial(selrpt);
  data.time   = data.time(selrpt);
  if isfield(data, 'offset'), data.offset = data.offset(selrpt); end

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

if ~isempty(selchan)
  for k = 1:numel(data.trial)
    data.trial{k} = data.trial{k}(selchan, :);
  end
  data.label = data.label(selchan);
end
