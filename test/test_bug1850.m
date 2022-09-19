function test_bug1850

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_neighbours ft_channelrepair
%
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1850

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));

cfg = [];
cfg.channel = {'all', '-MRT23', '-MLP57'};
data = ft_selectdata(cfg, data);

cfg = [];
cfg.method = 'template';
cfg.template = 'CTF275_neighb.mat';
neighbours = ft_prepare_neighbours(cfg);

% get the full list of 275 channel names
allchans = {neighbours.label};
missingchans = setdiff(allchans, data.label);

% repair the two channels that were removed
cfg = [];
cfg.missingchannel = missingchans;
cfg.neighbours = neighbours;
cfg.method = 'spline';
data_repaired = ft_channelrepair(cfg,data);
