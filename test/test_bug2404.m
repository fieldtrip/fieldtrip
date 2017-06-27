function test_bug2404

% WALLTIME 00:10:00
% MEM 1500mb

% TEST ft_channelrepair

data = [];
data.label = {'1', '2', '3'};
data.time{1} = (1:1000)/1000;
data.trial{1} = randn(3,1000);
data.trial{1}(2,:) = 0; % broken channel

elec.label = {'1', '2', '3', '4'};
elec.elecpos = randn(4,3);

grad.label = {'1', '2', '3', '4'};
grad.coilpos = randn(4,3);
grad.coilori = randn(4,3);

eegdata = data;
eegdata.elec = elec;

megdata = data;
megdata.grad = grad;

cfg = [];
cfg.badchannel = '2';
cfg.neighbours(1).label = '2';
cfg.neighbours(1).neighblabel = {'1', '3'};
fixedeeg = ft_channelrepair(cfg, eegdata);

cfg = [];
cfg.badchannel = '2';
cfg.neighbours(1).label = '2';
cfg.neighbours(1).neighblabel = {'1', '3'};
fixedmeg = ft_channelrepair(cfg, megdata);

assert(isfield(fixedeeg, 'elec'));
assert(isfield(fixedmeg, 'grad'));
assert(~all(fixedeeg.trial{1}(2,:)==0));
assert(~all(fixedmeg.trial{1}(2,:)==0));
