function test_bug2834

% MEM 1gb
% WALLTIME 00:20:00

% TEST ft_sourceanalysis

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% get some data
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg/freq_mtmconvol_powandcsd_ctf151.mat');
load(filename);

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_singleshell.mat');
load(filename);

% the issue seems to be due to the fact that ft_selectdata with cfg.channel
% does not automatically subselect the corresponding channelcmbs

cfg = [];
cfg.channel = {'MEG' '-MLO21' '-MLF21'};
freq2 = ft_selectdata(cfg, freq);

% this works from the point of view of ft_selectdata
tmp = repmat(ft_channelselection(cfg.channel, freq.label), [1 149]);
tmp2 = tmp';

cfg = [];
cfg.channel = {'MEG' '-MLO21' '-MLF21'};
cfg.channelcmb = [tmp(tril(ones(149),-1)>0) tmp2(tril(ones(149),-1)>0)];
freq3 = ft_selectdata(cfg, freq);

cfg = [];
cfg.vol = vol;
cfg.grid.resolution = 1;
cfg.method = 'dics';
cfg.frequency = 14;
cfg.latency   = 0.5;
source = ft_sourceanalysis(cfg, freq3);
