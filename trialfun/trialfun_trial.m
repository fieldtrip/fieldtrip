function [trl, event] = trialfun_trial(cfg);

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
  begsample = event(sel(i)).sample;
  endsample = begsample + event(sel(i)).duration - 1;
  offset    = event(sel(i)).offset;
  % add it to the list
  trl = [trl; [begsample endsample offset]];
end

