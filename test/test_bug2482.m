function test_bug2482

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_definetrial ft_preprocessing ft_trialfun_brainvision_segmented

%% define the trials

headerfile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2482/001 P11 M1 004_Seg_fs6.vhdr');
datafile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2482/001 P11 M1 004_Seg_fs6.eeg');

cfg = [];
cfg.trialfun = 'ft_trialfun_brainvision_segmented';
cfg.trigformat = 'DIN%d';

cfg.dataset = headerfile;
cfg1 = ft_definetrial(cfg);

cfg.dataset = datafile;
cfg2 = ft_definetrial(cfg);

if ~isequal(cfg1.trl, cfg2.trl)
  error('different trl matrix generated when specifying data file or headerfile as cfg.dataset');
end

%% read data

data = ft_preprocessing(cfg1);

%% test to see whether ERP is fine
% no need to run this every time, just a one time check (ES 3-mar-2014)
% cfg = [];
% cfg.channel = 'Cz';
% cfg.vartrllength = 2;
% tl = ft_timelockanalysis(cfg, data);
% 
% cfg = [];
% cfg.ylim = [-50 50];
% ft_singleplotER(cfg,tl);
