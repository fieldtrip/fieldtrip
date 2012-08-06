function [trl, event] = trialfun_trial(cfg)

% TRIALFUN_TRIAL creates a trial definition that corresponds to the
% events that are returned by FT_READ_EVENT with type='trial'
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

if isfield(cfg, 'event')
  % for BCI applications events should be specified in the cfg
  % to prevent reading the same events many times
  event = cfg.event;
else
  event = ft_read_event(cfg.dataset);
end

sel = find(strcmp({event.type}, 'trial'));
trl = zeros(length(sel),3);

smp = [event.sample];

for i=1:length(sel)
  % determine the begin, end and offset for each trial and add it to the Nx3 matrix
  begsample = event(sel(i)).sample;
  endsample = begsample + event(sel(i)).duration - 1;
  offset    = event(sel(i)).offset;

  tmpsel    = find(smp==smp(sel(i))-offset);
  tmpval    = zeros(1,0);
  for k=1:length(tmpsel)
    if ~isempty(event(tmpsel(k)).value)
      % this assumes per trial that the triggers occur in a standardised
      % order, otherwise the entries per column will have different
      % meanings
      tmpval = [tmpval event(tmpsel(k)).value];
    end
  end
  
  trl(i,1:(length(tmpval)+3))  = [begsample endsample offset tmpval];
end

