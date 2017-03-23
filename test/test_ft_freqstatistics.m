function test_ft_freqstatistics

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqstatistics, findcluster, clusterstat, ft_statistics_montecarlo

global ft_default;
ft_default.feedback = 'no';

% copyright, Roemer, bug 1201 (copyright? really? did I really put this in here? :P) - roevdmei

%%***********************************
% Expected support for freq_time without chan was removed by Robert around r9309, see
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1201#c20
%%***********************************
% For the case of "freq_time"
%
% % make fake dataset
% freq = cell(1,10);
% for idat = 1:10
%   freq{idat}.label = {'singlechan'};
%   freq{idat}.dimord = 'freq_time';
%   freq{idat}.powspctrm = rand(10,30);
%   freq{idat}.freq = 1:10;
%   freq{idat}.time = 0.1:0.1:3;
%   freq{idat}.cfg = [];
% end
% 
% % do stats - montecarlo
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05; 
% cfg.correctm    = 'cluster'; 
% cfg.clusterstatistic = 'maxsum';
% cfg.clusterthreshold = 'parametric';
% cfg.numrandomization = 500;
% cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
% cfg.ivar   = 1;
% cfg.uvar   = 2;
% stat = ft_freqstatistics(cfg,freq{:});
% 
% 
% % do stats - analytic
% cfg = [];
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05; 
% cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
% cfg.ivar   = 1;
% cfg.uvar   = 2;
% stat = ft_freqstatistics(cfg,freq{:});
% 
% 
% 
% % do stats - stats
% cfg = [];
% cfg.method      = 'stats';
% cfg.statistic   = 'ttest';
% cfg.alpha       = 0.05; 
% cfg.design = [ones(1,10) ];
% stat = ft_freqstatistics(cfg,freq{:});


%%***********************************
%%***********************************
%%***********************************
% And for the case of 'chan_freq_time' with one channel

% make fake dataset
freq = cell(1,10);
for idat = 1:10
  freq{idat}.label = {'singlechan'};
  freq{idat}.dimord = 'chan_freq_time';
  freq{idat}.powspctrm = rand(1,10,30);
  freq{idat}.freq = 1:10;
  freq{idat}.time = 0.1:0.1:3;
  freq{idat}.cfg = [];
end

% do stats - montecarlo
cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.correctm    = 'cluster'; 
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.numrandomization = 500;
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;
stat = ft_freqstatistics(cfg,freq{:});

% do stats - analytic
cfg = [];
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05; 
cfg.design = [ones(1,5) ones(1,5).*2; 1:5 1:5;];
cfg.ivar   = 1;
cfg.uvar   = 2;
stat = ft_freqstatistics(cfg,freq{:});



% do stats - stats
cfg = [];
cfg.method      = 'stats';
cfg.statistic   = 'ttest';
cfg.alpha       = 0.05; 
cfg.design = [ones(1,10) ];
stat = ft_freqstatistics(cfg,freq{:});



