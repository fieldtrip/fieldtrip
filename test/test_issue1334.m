function test_issue1334

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_selectdata ft_channelselection

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_issue1334.mat'));

%%

cfg = [];
cfg.channel={'megref'}; % reference data

dataRef = ft_preprocessing(cfg, data);
assert(numel(dataRef.label)==27);

%%
% see https://github.com/fieldtrip/fieldtrip/pull/1342#issuecomment-601265279

data = [];
data.label = {};

% this should return empty, but it gave an error
assert(isempty(ft_channelselection('MLO', data)));
