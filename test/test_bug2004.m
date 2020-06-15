function test_bug2004

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivityanalysis ft_connectivity_corr ft_connectivity_powcorr_ortho

%% test the functionality of ft_connectivityanalysis with respect to source level data (pcc)

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2004.mat');
load(filename);

cfg = [];
cfg.method = 'coh';
out = ft_connectivityanalysis(cfg, data);
%cfg.method = 'plv';
%out2 = ft_connectivityanalysis(cfg, data);
cfg.method = 'amplcorr';
out3 = ft_connectivityanalysis(cfg, data);
cfg.method = 'powcorr';
out4 = ft_connectivityanalysis(cfg, data);
cfg.method = 'powcorr_ortho';
out5 = ft_connectivityanalysis(cfg, data);

% grabbing the data from another bug.
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2468.mat');
load(filename);

% recompute source level pcc data
cfg                             = [];
cfg.frequency                   = fdata.freq;
cfg.headmodel                   = sourceVol;
cfg.sourcemodel                 = leadfieldModel;
cfg.sourcemodel.filter          = spatialFilters.avg.filter;
cfg.sourcemodel.label           = spatialFilters.avg.label;
cfg.keeptrials                  = 'no';
cfg.method                      = 'pcc';
cfg.(cfg.method).keepfilter     = 'yes';
cfg.(cfg.method).lambda         = '5%';
cfg.(cfg.method).fixedori       = 'no';
cfg.(cfg.method).keepmom        = 'yes';
sdatanew                        = ft_sourceanalysis(cfg, fdata);

% now also test the new 'functionality' for the other options of
% ft_sourcedescriptives
cfg            = [];
cfg.projectmom = 'yes';
cfg.keeptrials = 'yes';
cfg.fixedori   = 'over_trials';
sdnew          = ft_sourcedescriptives(cfg, sdatanew);

cfg = [];
cfg.method = 'coh';
out = ft_connectivityanalysis(cfg, sdnew);
%cfg.method = 'plv';
%out2 = ft_connectivityanalysis(cfg, data);
cfg.method = 'amplcorr';
out3 = ft_connectivityanalysis(cfg, sdnew);
cfg.method = 'powcorr';
out4 = ft_connectivityanalysis(cfg, sdnew);
cfg.method = 'powcorr_ortho';
out5 = ft_connectivityanalysis(cfg, sdnew);

