function test_bug1232

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_selectdata ft_datatype_source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

source = [];
source.pos = randn(100,3);
source.inside  = 1:2:100;
source.outside = 2:2:100;
for i=1:10
  source.trial(i).pow = nan(100,1); % it should be a column vector
  source.trial(i).pow(source.inside) = i*ones(50,1);
end


cfg = [];
cfg.trials = [2 4 6 8 10];
output = ft_selectdata(cfg, source);

assert(isequal(output.pow(1,:), cfg.trials), 'incorrect selection of trials in source structure');
assert(all(isnan(output.pow(2,:))),          'incorrect selection of trials in source structure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

source = [];
source.pos = randn(100,3);
source.inside  = 1:2:100;
source.outside = 2:2:100;
for i=1:10
  source.trial(i).pow = cell(100,1);
  for j=1:length(source.inside)
    source.trial(i).pow{source.inside(j)} = i;
  end
end


cfg = [];
cfg.trials = [2 4 6 8 10];
output = ft_selectdata(cfg, source);

assert(isequal(output.pow{1}, cfg.trials(:)), 'incorrect selection of trials in source structure');
assert(isempty(output.pow{2}),                'incorrect selection of trials in source structure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

source = [];
source.pos = randn(100,3);
source.inside  = 1:2:100;
source.outside = 2:2:100;
for i=1:10
  source.trial(i).csd = cell(100,1);
  for j=1:length(source.inside)
    source.trial(i).csd{source.inside(j)} = i*eye(3);
  end
end

cfg = [];
cfg.trials = [2 4 6 8 10];
output = ft_selectdata(cfg, source);

assert(isequal(squeeze(output.csd{1}(1,:,:)), 2*eye(3)), 'incorrect selection of trials in source structure');
assert(isempty(output.csd{2}),                           'incorrect selection of trials in source structure');

