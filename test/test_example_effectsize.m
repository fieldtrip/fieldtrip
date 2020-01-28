function test_example_effectsize

% WALLTIME 00:20:00
% MEM 4gb
% DEPENDENCY ft_statistics_analytic ft_statfun_cohensd


%%

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/effectsize/ERF_orig.mat'))


%%

cfg = [];
cfg.keepindividual = 'yes';
grandavgFIC = ft_timelockgrandaverage(cfg, allsubjFIC{:});
grandavgFC  = ft_timelockgrandaverage(cfg, allsubjFC{:});


%%

cfg = [];
cfg.channel = 'MLT12';
cfg.latency = [0.3 0.7];
cfg.avgoverchan = 'no';  % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
roiFIC = ft_selectdata(cfg, grandavgFIC);
roiFC  = ft_selectdata(cfg, grandavgFC);

x1 = roiFIC.individual(:)*1e15; % express the data in fT
x2 = roiFC.individual(:)*1e15;  % express the data in fT

n1 = length(x1);
n2 = length(x2);

figure; plot([x1 x2]', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
  'subj7', 'subj8', 'subj9', 'subj10'}, 'location', 'EastOutside');

%% non-paired, this is not optimal since it is a paired experimental design

pooled_sd = sqrt( ((n1-1)*std(x1)^2 + (n2-1)*std(x2)^2) / (n1+n2-1) );
cohensd = abs(mean(x1)-mean(x2)) / pooled_sd;

disp(cohensd)

% see https://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d

% Very small  0.01
% Small       0.20
% Medium      0.50
% Large       0.80
% Very large  1.20
% Huge        2.00


%%

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  ];
effect_roi_unpaired = ft_timelockstatistics(cfg, roiFIC, roiFC);

disp(effect_roi_unpaired)


%% paired

cohensd = mean(x1-x2) ./ std(x1-x2);
disp(cohensd)


%%

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_roi = ft_timelockstatistics(cfg, roiFIC, roiFC);

disp(effect_roi)


%% we can also do this for all channels and time points

cfg = [];
cfg.parameter = 'individual';
cfg.channel = 'MEG';
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_all = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);

% we can determine the channel and latency with the maximum effect
[m, ind] = max(effect_all.cohensd(:));
[i, j]   = ind2sub(size(effect_all.cohensd), ind);
fprintf('\nmaximum effect of %g on channel %s at latency %g\n\n', effect_all.cohensd(i,j), effect_all.label{i}, effect_all.time(j));

% or plot a distribution of the effect over all channels
cfg = [];
cfg.layout = 'CTF151_helmet.mat';
cfg.parameter = 'cohensd';
ft_multiplotER(cfg, effect_all);


%% we can also compute the effect for an average of the data in a pre-specified ROI (channels and time)

cfg = [];
cfg.channel = {'MLT12', 'MLT13', 'MLT23', 'MLT24', 'MLT32', 'MLT33', 'MLT41'};
cfg.latency = [0.35 0.55];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_avg = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);

disp(effect_avg)


%%

% define neighbouring channels
load ctf151_neighb.mat

cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';

cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'no';

cfg.alpha       = 0.05;
cfg.numrandomization = 1000;
cfg.neighbours  = neighbours; % defined as above
cfg.minnbchan   = 2; % minimal neighbouring channels
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.spmversion  = 'spm12';

cfg.ivar   = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar   = 2; % the 2nd row in cfg.design contains the subject number
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];

cfg.parameter   = 'avg';
inference = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});


%%

close all

figure; imagesc(inference.time, 1:151, -log10(inference.prob)); colorbar

% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = find(any(inference.mask,2));
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, grandavgFIC)
title('Nonparametric test with cluster-based multiple comparison correction')


%%

% % make a plot
% cfg = [];
% cfg.highlightsymbolseries = ['*','*','.','.','.'];
% cfg.layout = 'CTF151_helmet.mat';
% cfg.contournum = 0;
% cfg.markersymbol = '.';
% cfg.alpha = 0.05;
% cfg.parameter='stat';
% cfg.zlim = [-5 5];
% cfg.saveaspng = 'inference';
% ft_clusterplot(cfg, inference);


%%

close all

cfg = [];
cfg.channel = 'MEG';
cfg.latency = [0.3 0.7];
effect_with_mask = ft_selectdata(cfg, effect_all);
effect_with_mask.mask = inference.mask;


cfg = [];
cfg.layout = 'CTF151_helmet.mat';
cfg.parameter = 'cohensd';
cfg.maskparameter = 'mask';
ft_multiplotER(cfg, effect_with_mask);


%%

cfg = [];
cfg.channel     = 'MEG';
cfg.latency     = [0.300 0.700];
grandavgFIC_sel = ft_selectdata(cfg, grandavgFIC);
grandavgFC_sel  = ft_selectdata(cfg, grandavgFC);

x1 = nan(10,1);
x2 = nan(10,1);

for i=1:10
  sel3d = false(size(grandavgFIC_sel.individual));
  sel3d(i,:,:) = inference.mask;
  
  tmp = grandavgFIC_sel.individual(sel3d(:));
  x1(i) = mean(tmp);

  tmp = grandavgFC_sel.individual(sel3d(:));
  x2(i) = mean(tmp);
end


n1 = length(x1);
n2 = length(x2);

figure; plot([x1 x2]', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
  'subj7', 'subj8', 'subj9', 'subj10'}, 'location', 'EastOutside');
title('individual scores, averaged over cluster');

cohensd = mean(x1-x2) ./ std(x1-x2);
disp(cohensd)


%%

cfg = [];
cfg.channel = inference.label(any(inference.mask,2));
cfg.latency = [min(inference.time(any(inference.mask,1))) max(inference.time(any(inference.mask,1)))];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_rectangle = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);
disp(effect_rectangle)


