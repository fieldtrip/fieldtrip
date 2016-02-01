function test_bug2469

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_connectivityanalysis
% TEST ft_connectivity_powcorr_ortho

% this test function tests the functionality of Hipp's method on
% sensor-level data.

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2469.mat');
load(filename);

% this is only 1 channel data, will not work to begin with.
fdata.label{2,1} = strrep(fdata.label{1},'left','right');
fdata.fourierspctrm(:,2) = fdata.fourierspctrm(:,1).*exp(1i.*rand(8,1));

cfg        = [];
cfg.method = 'powcorr_ortho';
x          = ft_connectivityanalysis(cfg, fdata);
  



