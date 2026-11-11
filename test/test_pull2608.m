function test_pull2608

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_statistics_montecarlo clusterstat resampledesign ft_statfun_depsamplesT
% DATA no

% This test covers pull request 2608 (issue 2607): Monte Carlo p-values must
% count randomizations whose statistic is *at least as extreme* as the observed
% one (>=), not strictly more extreme (>). It checks:
%
%   1) with cfg.numrandomization='all' the identity permutation is a tie with
%      the observed statistic and must be counted, so the smallest possible
%      p-value is 1/Nperm (previously 0) and every p-value equals
%      #{T >= tobs}/Nperm computed independently;
%   2) the same for cfg.correctm='max';
%   3) cluster p-values are consistent with the returned randomization
%      distribution under the >= definition, including for the integer-valued
%      'maxsize' statistic where ties are frequent.

% ----------------------------------------------------------------------------
% 1) uncorrected p-values with all 2^10 sign flips of a paired design
% ----------------------------------------------------------------------------
rng(1,'twister');
nsubj  = 10;
design = [ones(1,nsubj) 2*ones(1,nsubj); 1:nsubj 1:nsubj];

cfg = [];
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.ivar             = 1;
cfg.uvar             = 2;
cfg.numrandomization = 'all';   % 1024 sign flips, includes the identity
cfg.correctm         = 'no';
cfg.tail             = 1;
cfg.feedback         = 'no';

% a strong effect: the observed t is the most extreme of the enumeration
dat = randn(1, 2*nsubj);
dat(1:nsubj) = dat(1:nsubj) + 3;
stat = ft_statistics_montecarlo(cfg, dat, design);
assert(abs(stat.prob - 1/1024) < 1e-12, ...
  'smallest possible p-value with all permutations should be 1/1024, got %g', stat.prob);

% a weaker effect: compare with an independent enumeration of the null distribution
dat2 = dat; dat2(1:nsubj) = dat2(1:nsubj) - 2.4;
stat2 = ft_statistics_montecarlo(cfg, dat2, design);
d    = dat2(1:nsubj) - dat2(nsubj+1:end);
tobs = sqrt(nsubj)*mean(d)/std(d);
T    = zeros(1, 2^nsubj);
for i = 1:2^nsubj
  s    = 1 - 2*(dec2bin(i-1,nsubj)=='1');
  ds   = d.*s;
  T(i) = sqrt(nsubj)*mean(ds)/std(ds);
end
pref = sum(T >= tobs)/2^nsubj;
assert(abs(stat2.prob - pref) < 1e-12, ...
  'p-value with all permutations should be #{T>=tobs}/Nperm = %g, got %g', pref, stat2.prob);

% ----------------------------------------------------------------------------
% 2) max-statistic correction with all permutations: p >= 1/Nperm everywhere
% ----------------------------------------------------------------------------
cfg.correctm = 'max';
dat3 = randn(5, 2*nsubj); dat3(:,1:nsubj) = dat3(:,1:nsubj) + 3;
stat3 = ft_statistics_montecarlo(cfg, dat3, design);
assert(all(stat3.prob(:) >= 1/1024 - 1e-12), ...
  'max-statistic p-values with all permutations must not be smaller than 1/Nperm');

% ----------------------------------------------------------------------------
% 3) cluster p-values agree with the >= definition on the returned distribution
% ----------------------------------------------------------------------------
T = 60; nsubj = 12; nrand = 300;
design = [ones(1,nsubj) 2*ones(1,nsubj); 1:nsubj 1:nsubj];
cfg = [];
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.ivar             = 1;
cfg.uvar             = 2;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clustertail      = 1;
cfg.tail             = 1;
cfg.numrandomization = nrand;
cfg.feedback         = 'no';
cfg.dim              = [1 T];
cfg.dimord           = 'chan_time';
cfg.channel          = {'c1'};
cfg.neighbours       = [];

for cs = {'maxsize', 'maxsum'}
  cfg.clusterstatistic = cs{1};
  x   = randn(T+10, 2*nsubj);
  x   = filter(ones(1,8)/8, 1, x);  % temporal smoothing so that clusters form
  dat = x(11:end,:);
  dat(20:30, 1:nsubj) = dat(20:30, 1:nsubj) + 0.6;
  stat = ft_statistics_montecarlo(cfg, dat, design);
  assert(isfield(stat, 'posclusters') && ~isempty(stat.posclusters), ...
    'expected at least one positive cluster for clusterstatistic=%s', cs{1});
  for j = 1:numel(stat.posclusters)
    s    = stat.posclusters(j).clusterstat;
    pref = (sum(stat.posdistribution >= s) + 1)/(nrand + 1);
    assert(abs(stat.posclusters(j).prob - pref) < 1e-12, ...
      '%s cluster %d: prob %g differs from >= definition %g', cs{1}, j, stat.posclusters(j).prob, pref);
  end
end
