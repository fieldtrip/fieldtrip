function [data] = selfromraw(data, varargin)

selrpt  = keyval('rpt',  varargin{:}); selectrpt  = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'rpt'));
selchan = keyval('chan', varargin{:}); selectchan = ~isempty(strmatch(varargin(cellfun(@ischar, varargin)), 'chan'));

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
