function test_tutorial_eventrelatedstatistics(datadir)

% MEM 2gb
% WALLTIME 00:20:00

% TEST test_tutorial_eventrelatedstatistics
% TEST ft_timelockanalysis ft_multiplotER ft_singleplotER ft_timelockstatistics
% TEST ft_topoplotER ft_clusterplot

% disable verbose output
global ft_default;
ft_default.feedback = 'no';

if nargin==0
  if ispc
    datadir = 'H:';
  else
    datadir = '/home';
  end
  
  load(fullfile(datadir, 'common', 'matlab', 'fieldtrip', 'data', 'ftp', 'tutorial', 'cluster_permutation_timelock', 'GA_ERF_orig.mat'));
else
  load(fullfile(datadir, 'GA_ERF_orig.mat'));
end


 %% plotting the grand-average and the subject-averages
 cfg = [];
cfg.showlabels  = 'yes';
cfg.layout    	= 'CTF151.lay';
figure; ft_multiplotER(cfg,GA_FC, GA_FIC)
 
cfg = [];
cfg.channel = 'MLT12';
figure; ft_singleplotER(cfg,GA_FC, GA_FIC)

%
figure; 
for iSub = 1:10
  subplot(3,4,iSub)
  plot(GA_FC.time,squeeze(GA_FC.individual(iSub,52,:)))
  hold on
  plot(GA_FIC.time,squeeze(GA_FIC.individual(iSub,52,:)),'r')
  title(strcat('subject ',num2str(iSub)))
  ylim([0 1.9e-13])
  xlim([-1 2])
end
subplot(3,4,11); 
text(0.5,0.5,'FC','color','b') ;text(0.5,0.3,'FIC','color','r')
axis off

%
chan = 52;
time = [0.3 0.7];
 
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
values_FIC = mean(GA_FIC.individual(:,chan,timesel_FIC),3);
 
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
values_FC = mean(GA_FC.individual(:,chan,timesel_FC),3);
 
% plot to see the effect in each subject
M = [values_FC,values_FIC];
figure; plot(M','o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');
      
%% T-test with Matlab function

%dependent samples ttest
FCminFIC = values_FC - values_FIC;
[h,p,ci,stats] = ttest(FCminFIC, 0, 0.05) % H0: mean = 0, alpha 0.05

%% T-test with FieldTrip function
cfg = [];
cfg.channel     = 'MLT12';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'individual';
cfg.method      = 'analytic';
cfg.statistic   = 'depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,GA_FC,GA_FIC)

%% Multiple comparisons

%loop over channels
time = [0.3 0.7];
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
clear h p
for iChan = 1:151
  FICminFC = mean(GA_FIC.individual(:,iChan,timesel_FIC),3) - mean(GA_FC.individual(:,iChan,timesel_FC),3);
  [h(iChan), p(iChan)] = ttest(FICminFC, 0, 0.05 ); % test each channel separately
end
 
% plot uncorrected "significant" channels
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(h);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_FC)
title('significant without multiple comparison correction')


%% with Bonferoni correction for multiple comparisons
for iChan = 1:151
  FICminFC = mean(GA_FIC.individual(:,iChan,timesel_FIC),3) - mean(GA_FC.individual(:,iChan,timesel_FC),3);
  [h(iChan), p(iChan)] = ttest(FICminFC, 0, 0.05/151); % Bonferoni correction for 151 channels
end



%% Bonferoni correction with FieldTrip function
cfg = [];
cfg.channel     = 'MEG'; %now all channels
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'individual';
cfg.method      = 'analytic';
cfg.statistic   = 'depsamplesT'
cfg.alpha       = 0.05;
cfg.correctm    = 'bonferoni';
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,GA_FIC,GA_FC)

%% Permutation test based on t statistics
fg = [];
cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'individual';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT'
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,GA_FIC,GA_FC)
 
% make the plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_FC)
title('Nonparametric: significant without multiple comparison correction')

%%  Permutation test based on cluster statistics
cfg = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.feedback    = 'yes'; % show a neighbour plot 
neighbours      = ft_prepare_neighbours(cfg, GA_FC); % define neighbouring channels

cfg = [];
cfg.neighbours  = neighbours; % define neighbouring channels
cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'individual';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT'
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,GA_FIC,GA_FC)
 
% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_FC)
title('Nonparametric: significant with cluster multiple comparison correction')

% different timewindow
cfg = [];
cfg.neighbours  = neighbours; % define neighbouring channels
cfg.channel     = 'MEG';
cfg.latency     = [-0.5 2];
cfg.avgovertime = 'no';
cfg.parameter   = 'individual';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
cfg.minnbchan        = 2; % minimal neighbouring channels
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 

failed = true;
for i=1:100 % this can fail sometime, like every second iteration or so, 100 is a *really* conservative number here
  stat = ft_timelockstatistics(cfg,GA_FIC,GA_FC);
  
  % make a plot
  pcfg = [];
  pcfg.highlightsymbolseries = ['*','*','.','.','.'];
  pcfg.layout = 'CTF151.lay';
  pcfg.contournum = 0;
  pcfg.markersymbol = '.';
  pcfg.parameter = 'stat';
  pcfg.alpha = 0.05;
  pcfg.zlim = [-5 5];
  try
    ft_clusterplot(pcfg,stat);
    failed = false;
    break;
  end
end


if failed
  error('ft_clusterplot fails cause p was too high');
end
