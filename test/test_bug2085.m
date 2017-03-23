function test_bug2085

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_senstype ft_senslabel

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%% create a volume conductor
vol = [];
vol.r = 10;
vol.o = [0 0 0];
vol.unit = 'cm';
vol = ft_datatype_headmodel(vol);

%% create a set of sensors
[pnt, tri] = icosahedron162;
pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);
grad.pnt = pnt(sel,:) .* 1.2;
grad.ori = pnt(sel,:);
grad.tra = eye(length(sel));
for i=1:length(sel)
  grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
  grad.label{i} = sprintf('magnetometer%d', i);
end
grad.unit = 'cm';

grad1 = ft_datatype_sens(grad);
grad2 = ft_datatype_sens(grad);

grad2.type = 'magnetometer'; % this makes ft_senstype much faster

%% determine the time to compute some leadfields

cfg = [];
cfg.vol = vol;
cfg.grid.resolution = 4;
cfg.channel = 'all';

cfg.grad = grad1;
grid1 = ft_prepare_leadfield(cfg); % don't count the first one
tic
for i=1:10
  grid1 = ft_prepare_leadfield(cfg); % iterate a few times for better time estimate
end
time1 = toc;

cfg.grad = grad2;
grid2 = ft_prepare_leadfield(cfg); % don't count the first one
tic
for i=1:10
  grid2 = ft_prepare_leadfield(cfg); % iterate a few times for better time estimate
end
time2 = toc;

%% compare the time that it took, this fraction 1.5 is rather arbitrary

if (time1/time2)>1.5
  error('the leadfield computation without grad.type takes too long (%d seconds compared to %d seconds)', time1, time2);
end
