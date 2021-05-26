function test_bug2005

% WALLTIME 00:20:00
% MEM 4gb
% DEPENDENCY ft_sourceanalysis

fname = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2005.mat');
load(fname);

foi = 5;

cfg              = [];
cfg.foi          = foi;
cfg.trials       = 'all';
cfg.keeptrials   = 'yes';
cfg.keeptapers   = 'yes';
cfg.output       = 'fourier';
cfg.channel      = 'MEG';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 2;

f_data = ft_freqanalysis(cfg, data);
% f = frequency

% JM hack
cf_data = f_data;
cf_data = ft_checkdata(cf_data, 'cmbstyle', 'fullfast'); % efficiently compute crsspctrm
cf_data = ft_checkdata(cf_data, 'cmbstyle', 'sparse'); % remove redundancy in channelcomb
cf_data = ft_checkdata(cf_data, 'cmbstyle', 'sparsewithpow'); % restore power data
% data is now ready to be fed into sourceanalysis

% Source analysis
% Compute common spatial filter

cfg                             = [];
cfg.method                      = 'dics';
cfg.sourcemodel                 = sourceModelGrid;
cfg.headmodel                   = sourceModelVol;
cfg.frequency                   = foi;
cfg.keeptrials                  = 'yes';
cfg.(cfg.method).keepfilter     = 'yes';
cfg.(cfg.method).fixedori       = 'yes';
cfg.(cfg.method).realfilter     = 'yes';
cfg.(cfg.method).lambda         = '5%';

sf_data                         = ft_sourceanalysis(cfg, cf_data);
% s = source


cfg                             = [];
cfg.method                      = 'pcc';

cfg.sourcemodel                 = sourceModelGrid;
cfg.sourcemodel.filter          = sf_data.avg.filter;
cfg.headmodel                   = sourceModelVol;
cfg.frequency                   = foi;

cfg.keeptrials                  = 'yes';
cfg.(cfg.method).keepmom        = 'yes';

% error is produced here:
sf_data2                         = ft_sourceanalysis(cfg, f_data);
