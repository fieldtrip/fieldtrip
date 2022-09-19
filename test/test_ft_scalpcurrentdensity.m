function test_ft_scalpcurrentdensity

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_scalpcurrentdensity

fs = 500;
nchan = 32;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));

% construct a set of electrodes randomly distributed over the upper hemisphere
data.elec.label = data.label;
data.elec.unit = 'cm';
data.elec.tra = eye(nchan);
data.elec.elecpos = randn(nchan,3);
data.elec.elecpos(:,3) = abs(data.elec.elecpos(:,3));
for i=1:nchan
  data.elec.elecpos(i,:) = 10*data.elec.elecpos(i,:)/norm(data.elec.elecpos(i,:));
end
data.elec.chanpos = data.elec.elecpos;

cfg = [];
dataout = ft_scalpcurrentdensity(cfg, data);
