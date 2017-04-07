function test_bug2372

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_timelockgrandaverage ft_timelockanalysis

global ft_default;
ft_default.feedback = 'no';

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFC_LP.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFIC_LP.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataIC_LP.mat'));

%% creating timelock structures with avg and trial fields
cfg = [];
cfg.keeptrials = 'yes';
tlk1 = ft_timelockanalysis(cfg,dataFC_LP);
tlk2 = ft_timelockanalysis(cfg,dataFIC_LP);
tlk3 = ft_timelockanalysis(cfg,dataIC_LP);

%% checking for a general ft_timelockgrandaverage
%here ft_timelockgrandaverage should discard the trial field
cfg = [];
cfg.parameter = 'avg';
gavg = ft_timelockgrandaverage(cfg,tlk1,tlk2,tlk3);
  
%% test the cfg.parameter = 'trial' expected error
try
  cfg.parameter = 'trial';
  gavg = ft_timelockgrandaverage(cfg,tlk1,tlk2,tlk3);
  success=true;
catch
  success=false;
end
if success
  error('this should fail with cfg.parameter = trial');
end
