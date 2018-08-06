function test_bug843

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_topoplotTFR

% it has been reported that the linearly indexed connectivity metrics don't
% behave robustly in combination with a specified refchannel; also
% directionality needs to be documented

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275'));

% constrain to MEG channels
cfg = [];
cfg.channel = 'MEG';
data = ft_preprocessing(cfg, data);

% create a connectivity data structure
freq.labelcmb  = [data.label repmat(data.label(1),[numel(data.label) 1])];
freq.labelcmb  = freq.labelcmb(2:end,:);
freq.cohspctrm = rand(274,20);
freq.freq      = 1:20;
freq.dimord    = 'chancmb_freq';

cfg            = [];
cfg.layout     = 'CTF275.lay';
cfg.refchannel = 'MLC11';
cfg.parameter  = 'cohspctrm';
cfg.directionality = 'inflow';
figure;ft_topoplotTFR(cfg, freq);

%Crash, correct behavior
%cfg.directionality = 'outflow';
%ft_topoplotTFR(cfg, freq);

freq.labelcmb  = freq.labelcmb(:,[2 1]);
%Crash, correct behavior
%cfg.directionality = 'inflow';
%ft_topoplotTFR(cfg, freq);
cfg.directionality = 'outflow';
figure;ft_topoplotTFR(cfg, freq);

freq.labelcmb(1:140,:) = freq.labelcmb(1:140,[2 1]);
cfg.directionality = '';
figure;ft_topoplotTFR(cfg,freq);
cfg.directionality = 'inflow';
figure;ft_topoplotTFR(cfg,freq);
cfg.directionality = 'outflow';
figure;ft_topoplotTFR(cfg,freq);


