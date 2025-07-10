function test_ft_freqstatistics

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqstatistics findcluster clusterstat ft_statistics_montecarlo
% DATA no

if ~ft_hastoolbox('stats')
  warning('the STATS toolbox is not available, skipping this test');
  return
end

% copyright, Roemer, bug 1201 (copyright? really? did I really put this in here? :P) - roevdmei

%%***********************************
% Expected support for freq_time without chan was removed by Robert around r9309, see
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1201#c20
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

%%
% do stats, using a random seed, check whether the same results occur
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
stat1 = ft_freqstatistics(cfg,freq{:});
stat2 = ft_freqstatistics(cfg,freq{:});

cfg.randomseed = 42;
stat3 = ft_freqstatistics(cfg,freq{:});
stat4 = ft_freqstatistics(cfg,freq{:});

assert(~isequal(stat1.prob, stat2.prob));
assert(~isequal(stat1.prob, stat3.prob));
assert(~isequal(stat1.prob, stat4.prob));
assert( isequal(stat3.prob, stat4.prob));

%%
% call ft_statistics_montecarlo directly, this is not recommended, and
% should throw a warning. Moreover,the code below verifies that setting the
% randomseed for this function is not going to do the trick

% we need to create a data matrix by hand
for k = 1:numel(freq)
  dat(:,k) = freq{k}.powspctrm(:);
end
cfg.dim = [1 10 30]; % needed for ft_statistics_montecarlo

stat5 = ft_statistics_montecarlo(cfg, dat, cfg.design);
stat6 = ft_statistics_montecarlo(cfg, dat, cfg.design);

rng(42);
stat7 = ft_statistics_montecarlo(cfg, dat, cfg.design);
rng(42);
stat8 = ft_statistics_montecarlo(cfg, dat, cfg.design);

assert(~isequal(stat5.prob, stat6.prob));
assert(~isequal(stat5.prob, stat7.prob));
assert(~isequal(stat5.prob, stat8.prob));
assert( isequal(stat7.prob, stat8.prob));
assert( isequal(stat3.prob(:), stat8.prob)); % this should be the same
