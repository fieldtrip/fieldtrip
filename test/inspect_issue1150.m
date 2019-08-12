function inspect_issue1150

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_rejectvisual

data = [];
data.label = {'Cz'};
data.fsample = 1000;
for i=1:100
  data.time{i} = (1:1000) * 1/1000;
  data.trial{i}(1,:) = randn(1,1000);
end

% please click on all metric options, notably on the "min"
cfg = [];
clean = ft_rejectvisual(cfg, data);
