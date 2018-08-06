function test_bug2874

% WALLTIME 00:10:100
% MEM 1000mb

% TEST ft_sourcegrandaverage

pos = randn(19344, 3);
nsubj = 3;

clear source
for i=1:nsubj
  source{i}.pos = pos; % needs to be the same over subjects
  source{i}.dim = [24 31 26];
  source{i}.unit = 'cm';
  source{i}.inside = ones(19344, 1);
  source{i}.params = 'something'; % don't know what this is
  source{i}.initial = randn(4,4); % don't know what this is
  source{i}.cfg = [];
  source{i}.seed_l     = randn(19344, 1);
  source{i}.avg.seed_r = randn(19344, 1);
end

%%

cfg = [];
cfg.parameter = 'seed_l';
grandavg = ft_sourcegrandaverage(cfg, source{:});

assert(isfield(grandavg, 'seed_l'));

%%

cfg = [];
cfg.parameter = 'avg.seed_r';
grandavg = ft_sourcegrandaverage(cfg, source{:});

assert(isfield(grandavg, 'seed_r')); % in the output it is not in the avg substructure any more
