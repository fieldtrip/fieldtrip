function test_bug1448

% MEM 1500mb
% WALLTIME 00:10:00


% this function tests whether the mask is kept inside the call to singleplotTFR
% load data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg/freq_mtmconvol_trl_ctf275.mat'));

% create mask-field
freqsize = size(freq.powspctrm);
freq.maskfieldname = zeros(freqsize(2:4));
freq.maskfieldname(:,1:round(freqsize(3)/2),:) = 1;

% plot using singleplotTFR
cfg = [];
cfg.channel = 1;
cfg.maskparameter = 'maskfieldname';
ft_singleplotTFR(cfg,freq)
