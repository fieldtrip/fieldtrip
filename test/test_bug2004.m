function test_bug2004

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2004 
% TEST ft_connectivityanalysis
% TEST ft_connectivity_corr
% TEST ft_connectivity_powcorr_ortho

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

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
cfg.vol                         = sourceVol;
cfg.grid                        = leadfieldModel;
cfg.grid.filter                 = spatialFilters.avg.filter;
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

