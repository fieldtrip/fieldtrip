function test_ft_multiplotTFR

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_multiplotTFR 
% TEST ft_freqanalysis ft_multiplotTFR ft_prepare_layout

% writeflag determines whether the output should be saved to disk
% version determines the output directory

% dataset taken from ref_datasets, using the tutorial one for testing
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));

% get freqs
cfg            = [];
cfg.channel    = {'MEG'};
cfg.method     = 'mtmconvol';
cfg.output     = 'pow';
cfg.keeptrials = 'no';
cfg.foi        = 2:2:30;
cfg.taper      = 'hanning';
cfg.t_ftimwin  = ones(1,numel(cfg.foi)).*0.5;
cfg.toi        = (250:50:750)./1000;
freq = ft_freqanalysis(cfg,data);


% plot them
cfg = [];
ft_multiplotTFR(cfg,freq)

