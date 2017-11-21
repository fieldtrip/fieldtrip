function test_bug3370

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_topoplotER test_bug3370

% this test function requires user interaction to judge whether it runs
% through well, so it is of limited use in a automatic test batch setting
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

for k = 1:numel(data.trial)
  data.trial{k}(end+(1:2),:) = randn(2,numel(data.time{k}));
end
data.label = [data.label;{'ref1';'ref2'}];

cfg = [];
cfg.trials = 1:5;
data1 = ft_selectdata(cfg, data);
cfg.trials = 6:10;
data2 = ft_selectdata(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper  = 'hanning';
cfg.channel = {'MEG';'ref1';'ref2'};
cfg.foilim = [0 20];
freq1 = ft_freqanalysis(cfg, data1);
freq2 = ft_freqanalysis(cfg, data2);

cfg = [];
cfg.method = 'coh';
cfg.channelcmb = [ft_channelcombination({'all' 'ref1'},freq1.label);ft_channelcombination({'all' 'ref2'},freq1.label)];
coh1 = ft_connectivityanalysis(cfg, freq1);
coh2 = ft_connectivityanalysis(cfg, freq2);

cfg = [];
cfg.layout = 'CTF151_helmet.mat';
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'ref1';
figure;hold on;
subplot(2,1,1);ft_topoplotER(cfg, coh1);
cfg.refchannel = 'ref2';
subplot(2,1,2);ft_topoplotER(cfg, coh2);

keyboard



