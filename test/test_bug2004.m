function test_bug2004

% MEM 1gb
% TEST test_bug2004 
% TEST ft_connectivityanalysis
% TEST ft_connectivity_corr
% TEST ft_connectivity_powcorr_ortho

%% test the functionality of ft_connectivityanalysis with respect to source level data (pcc)

load ('/home/common/matlab/fieldtrip/data/test/bug2004.mat');

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
