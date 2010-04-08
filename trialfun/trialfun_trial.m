function [trl, event] = trialfun_trial(cfg)

% TRIALFUN_TRIAL creates a trial definition that corresponds to the
% events that are returned by READ_EVENT with type='trial'

if isfield(cfg, 'event')
  % for BCI applications events should be specified in the cfg
  % to prevent reading the same events many times
  event = cfg.event;
else
  event = read_event(cfg.dataset);
end

trl = [];
sel = find(strcmp({event.type}, 'trial'));

for i=1:length(sel)
  % determine the begin, end and offset for each trial and add it to the Nx3 matrix
  begsample = event(sel(i)).sample;
  endsample = begsample + event(sel(i)).duration - 1;
  offset    = event(sel(i)).offset;
  trl       = [trl; [begsample endsample offset]];
end

