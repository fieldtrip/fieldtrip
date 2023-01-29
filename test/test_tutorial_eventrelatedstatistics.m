function test_tutorial_eventrelatedstatistics

% MEM 4gb
% WALLTIME 00:20:00
% DEPENDENCY ft_timelockanalysis ft_multiplotER ft_singleplotER ft_timelockstatistics ft_topoplotER ft_clusterplot

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedstatistics/ERF_orig.mat'));

%% calculating the grand-average

% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grandavgFIC   = ft_timelockgrandaverage(cfg, allsubjFIC{:});
grandavgFC    = ft_timelockgrandaverage(cfg, allsubjFC{:});
% "{:}" means to use data from all elements of the variable


%% plotting the grand-averagecfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = 'CTF151_helmet.mat';
ft_multiplotER(cfg, grandavgFIC, grandavgFC)

cfg = [];
cfg.channel = 'MLT12';
ft_singleplotER(cfg, grandavgFIC, grandavgFC)


%%

time = [0.3 0.7];
% Scaling of the vertical axis for the plots below
ymax = 1.9e-13;
figure;
for isub = 1:10
  subplot(3,4,isub)
  % use the rectangle to indicate the time range used later
  rectangle('Position',[time(1) 0 (time(2)-time(1)) ymax],'FaceColor',[0.7 0.7 0.7]);
  hold on;
  % plot the lines in front of the rectangle
  plot(allsubjFIC{isub}.time,allsubjFIC{isub}.avg(52,:), 'b');
  plot(allsubjFC{isub}.time,allsubjFC{isub}.avg(52,:), 'r');
  title(strcat('subject ',num2str(isub)))
  ylim([0 1.9e-13])
  xlim([-1 2])
end
subplot(3,4,11);
text(0.5,0.5,'FIC','color','b') ;text(0.5,0.3,'FC','color','r')
axis off


%% plotting single subject averages

chan = 52;
time = [0.3 0.7];

% find the time points for the effect of interest in the grand average data
timesel_FIC = find(grandavgFIC.time >= time(1) & grandavgFIC.time <= time(2));
timesel_FC  = find(grandavgFC.time >= time(1) & grandavgFC.time <= time(2));

% select the individual subject data from the time points and calculate the mean
for isub = 1:10
  valuesFIC(isub) = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));
  valuesFC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
end

% plot to see the effect in each subject
M = [valuesFIC(:) valuesFC(:)];
figure; plot(M', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
  'subj7', 'subj8', 'subj9', 'subj10'}, 'location', 'EastOutside');


%% T-test with MATLAB function

FICminFC = valuesFIC - valuesFC;
[h,p,ci,stats] = ttest(FICminFC, 0, 0.05) % H0: mean = 0, alpha 0.05


%% T-test with FieldTrip function

% define the parameters for the statistical comparison
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

stat = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});   % don't forget the {:}!


%% Multiple comparisons

% loop over channels
time = [0.3 0.7];
timesel_FIC = find(grandavgFIC.time >= time(1) & grandavgFIC.time <= time(2));
timesel_FC  = find(grandavgFC.time >= time(1) & grandavgFC.time <= time(2));
clear h p

FICminFC = zeros(1,10);

for iChan = 1:151
  for isub = 1:10
    FICminFC(isub) = ...
      mean(allsubjFIC{isub}.avg(iChan,timesel_FIC)) - ...
      mean(allsubjFC{isub}.avg(iChan,timesel_FC));
  end
  
  [h(iChan), p(iChan)] = ttest(FICminFC, 0, 0.05 ); % test each channel separately
end

% plot uncorrected "significant" channels
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = find(h);
cfg.comment   = 'no';
ft_topoplotER(cfg, grandavgFIC)
title('significant without multiple comparison correction')


%% with Bonferoni correction for multiple comparisons

% with Bonferroni correction for multiple comparisons
FICminFC = zeros(1,10);

for iChan = 1:151
  for isub = 1:10
    FICminFC(isub) = ...
      mean(allsubjFIC{isub}.avg(iChan,timesel_FIC)) - ...
      mean(allsubjFC{isub}.avg(iChan,timesel_FC));
  end
  
  [h(iChan), p(iChan)] = ttest(FICminFC, 0, 0.05/151); % test each channel separately
end


%% Bonferoni correction with FieldTrip function

cfg = [];
cfg.channel     = 'MEG'; % now all channels
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'bonferroni';

Nsub = 10;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});


%% NON-PARAMETRIC Permutation test based on t-statistics

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

stat = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

% make the plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
ft_topoplotER(cfg, grandavgFIC)
title('Nonparametric: significant without multiple comparison correction')


%%  NON-PARAMETRIC: Permutation test based on cluster statistics

cfg = [];
cfg.method      = 'template';                         % try 'distance' as well
cfg.template    = 'ctf151_neighb.mat';                % specify type of template
cfg.layout      = 'CTF151_helmet.mat';                % specify layout of sensors
cfg.feedback    = 'yes';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, grandavgFIC); % define neighbouring channels

% note that the layout and template fields have to be entered because at the earlier stage
% when ft_timelockgrandaverage is called the field 'grad' is removed. It is this field that
% holds information about the (layout of the) sensors.

cfg = [];
cfg.channel     = 'MEG';
cfg.neighbours  = neighbours; % defined as above
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

stat = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF151_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
ft_topoplotER(cfg, grandavgFIC)
title('Nonparametric: significant with cluster-based multiple comparison correction')


%% longer time window for NON-PARAMETRIC: permutation test based on cluster statistics

cfg = [];
cfg.channel     = 'MEG';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [-0.25 1];
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
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

stat = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

% make a plot
cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout = 'CTF151_helmet.mat';
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
cfg.zlim = [-5 5];
cfg.toi  = [-0.1:0.1:0.9];
ft_clusterplot(cfg, stat);

