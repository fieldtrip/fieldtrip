function test_bug2834

% MEM 1gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceanalysis
% DATA private

% get some data
filename = dccnpath('/project/3031000.02/test/latest/freq/meg/freq_mtmconvol_powandcsd_ctf151.mat');
load(filename);

filename = dccnpath('/project/3031000.02/test/latest/vol/Subject01vol_singleshell.mat');
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
cfg.headmodel = vol;
cfg.sourcemodel.resolution = 1;
cfg.method = 'dics';
cfg.frequency = 14;
cfg.latency   = 0.5;
source = ft_sourceanalysis(cfg, freq3);
