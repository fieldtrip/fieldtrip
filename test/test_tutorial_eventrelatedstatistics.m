function test_tutorial_eventrelatedstatistics(datadir)

% MEM 2500mb
% WALLTIME 00:20:00

% TEST ft_timelockanalysis ft_multiplotER ft_singleplotER ft_timelockstatistics
% TEST ft_topoplotER ft_clusterplot

global ft_default;
ft_default.feedback = 'no';

if nargin==0
  if ispc
    datadir = 'H:';
  else
    datadir = '/home';
  end
  
  load(fullfile(datadir, 'common', 'matlab', 'fieldtrip', 'data', 'ftp', 'tutorial', 'eventrelatedstatistics', 'ERF_orig.mat'));
else
  load(fullfile(datadir, 'ERF_orig.mat'));
end

% Tutorial update on 3 Feb: no longer give grandaverage inarg to
% ft_XX_statistics, instead introduce cell-array inarg. 
%% calculating the grand-average

% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_FC         = ft_timelockgrandaverage(cfg,allsubjFC{:});  
GA_FIC        = ft_timelockgrandaverage(cfg,allsubjFIC{:});

%% plotting the grand-average
cfg = [];
cfg.showlabels  = 'yes';
cfg.layout    	= 'CTF151.lay';
figure; ft_multiplotER(cfg,GA_FC, GA_FIC)
 
cfg = [];
cfg.channel = 'MLT12';
figure; ft_singleplotER(cfg,GA_FC, GA_FIC)


time = [0.3 0.7];
% Scaling of the vertical axis for the plots below
ymax = 1.9e-13; 
figure; 
for iSub = 1:10
    subplot(3,4,iSub)
    plot(allsubjFC{iSub}.time,allsubjFC{iSub}.avg(52,:));
    hold on
    % use the rectangle to indicate the time range used later 
    rectangle('Position',[time(1) 0 (time(2)-time(1)) ymax],'FaceColor',[0.7 0.7 0.7]);
    % repeat the plotting to get the line in front of the rectangle
    plot(allsubjFIC{iSub}.time,allsubjFIC{iSub}.avg(52,:),'r');
    title(strcat('subject ',num2str(iSub)))
    ylim([0 1.9e-13])
    xlim([-1 2])
end
subplot(3,4,11); 
text(0.5,0.5,'FC','color','b') ;text(0.5,0.3,'FIC','color','r')
axis off

%% plotting single subject averages
chan = 52;
time = [0.3 0.7];
 
% find the time points for the effect of interest in the grand average data
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
 
% select the individual subject data from the time points and calculate the mean
for isub = 1:10
    values_FIC(isub)  = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));   
    values_FC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
end
 
% plot to see the effect in each subject
M = [values_FC',values_FIC'];
figure; plot(M','o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');
      
%% T-test with MATLAB function
% dependent samples ttest
FCminFIC = values_FC - values_FIC;
[h,p,ci,stats] = ttest_wrapper(FCminFIC, 0, 0.05) % H0: mean = 0, alpha 0.05;

%% T-test with FieldTrip function
cfg = [];
cfg.channel     = 'MLT12';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,allsubjFC{:},allsubjFIC{:});

%% Multiple comparisons

%loop over channels
time = [0.3 0.7];
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
clear h p
 
FICminFC = zeros(1,10);
 
for iChan = 1:151
    for isub = 1:10
        FICminFC(isub) = mean(allsubjFIC{isub}.avg(iChan,timesel_FIC)) - mean(allsubjFC{isub}.avg(iChan,timesel_FC));
        [h(iChan), p(iChan)] = ttest_wrapper(FICminFC, 0, 0.05 ); % test each channel separately
    end
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
FICminFC = zeros(1,10);
 
for iChan = 1:151
    for isub = 1:10
        FICminFC(isub) = mean(allsubjFIC{isub}.avg(iChan,timesel_FIC)) - mean(allsubjFC{isub}.avg(iChan,timesel_FC));
        [h(iChan), p(iChan)] = ttest_wrapper(FICminFC, 0, 0.05/151); % test each channel separately
    end
end


%% Bonferoni correction with FieldTrip function
cfg = [];
cfg.channel     = 'MEG'; %now all channels
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'bonferoni';
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
 
stat = ft_timelockstatistics(cfg,allsubjFC{:},allsubjFIC{:});


%% NON-PARAMETRIC Permutation test based on t statistics
cfg = [];
cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,allsubjFIC{:},allsubjFC{:});
 
% make the plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_FC)
title('Nonparametric: significant without multiple comparison correction')

%%  NON-PARAMETRIC: Permutation test based on cluster statistics
cfg = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.template    = 'ctf151_neighb.mat';               % specify type of template
cfg.layout      = 'CTF151.lay';                      % specify layout of sensors*
cfg.feedback    = 'yes';                             % show a neighbour plot 
neighbours      = ft_prepare_neighbours(cfg, GA_FC); % define neighbouring channels
 

cfg = [];
cfg.neighbours  = neighbours; % define neighbouring channels
cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,allsubjFIC{:},allsubjFC{:});
 
% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_FC)
title('Nonparametric: significant with cluster multiple comparison correction')

%% longer time window for NON-PARAMETRIC: permutation test based on cluster statistics

cfg = [];
cfg.channel     = 'MEG';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [0.1 0.5];
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
%cfg.tail        = 1;          % options are -1, 1 and 0 (default = 0, for 2-tailed)
cfg.numrandomization = 1000;
cfg.minnbchan        = 2; % minimal neighbouring channels
 
Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,allsubjFIC{:},allsubjFC{:});

failed = true;
for i=1:1000 % this can fail sometime, like every second iteration or so, 100 is a *really* conservative number here
  stat = ft_timelockstatistics(cfg,allsubjFIC{:},allsubjFC{:});
  
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


function [h,p,ci,stats]=ttest_wrapper(x,y,alpha)
% helper functions for ttest
% - old Matlab, with syntax:                   ttest(x,y,alpha,tail,dim)
% - new Matlab and GNU Octave, with syntax:    ttest(x,y,'alpha',alpha,...)

    if nargin('ttest')>0
        [h,p,ci,stats]=ttest(x,y,alpha);
    else
        [h,p,ci,stats]=ttest(x,y,'alpha',alpha);
    end
