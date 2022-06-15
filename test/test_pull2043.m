function test_pull2043

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY data2bids

ctfdataset = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043/emptyroom_Noise_20220613_01.ds');
neuromagdataset = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043/rest_ec_mc_avgtrans_tsss_corr95-raw.fif');

%%

bidsroot = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043'), 'bids');

rmdir(bidsroot, 's');

%%

cfg = [];
cfg.bidsroot = bidsroot;
cfg.dataset = ctfdataset;
cfg.datatype = 'meg';
cfg.method = 'copy';
cfg.sub = 'ctf';
cfg.ses = 'meg';
cfg.task = 'emptyroom';

cfg.participants.age = nan;
cfg.participants.sex = nan;

cfg.scans.acq_time = datetime('yesterday');

data2bids(cfg)

%%


cfg = [];
cfg.bidsroot = bidsroot;
cfg.dataset = neuromagdataset;
cfg.datatype = 'meg';
cfg.method = 'copy';
cfg.sub = 'neuromag';
cfg.ses = 'meg';
cfg.task = 'emptyroom';

cfg.participants.age = nan;
cfg.participants.sex = nan;

cfg.scans.acq_time = datetime('today');

data2bids(cfg)