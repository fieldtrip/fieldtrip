function test_ft_recodeevent

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_recodeevent

event = [];
trl = [];
for i = 1:100
  event(i).sample = i;
  event(i).offset = [];
  event(i).duration = 1;
  event(i).value = 1;
  trl(i,:) = [i,i+0.5,0.1];
end

cfg = [];
cfg.eventtype = [];
cfg.eventvalue = [];
eventout = ft_recodeevent(cfg,event,trl);
