function test_issue1334

% MEM 2gb
% WALLTIME 00:10:00

% DEPENDENCY ft_selectdata ft_channelselection

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_issue1334.mat'));

cfg = [];
cfg.channel={'megref'}; % reference data

dataRef = ft_preprocessing(cfg, data);
assert(numel(dataRef.label)==27);
