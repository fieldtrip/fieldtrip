function test_ft_math

% TEST test_ft_math
% TEST ft_math

% create some data
timelock1.label  = {'chan1'; 'chan2'};
timelock1.time   = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg    = ones(2,5);

timelock2 = timelock1;
timelock2.avg  = ones(2,5)*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with one input

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'log10';
tmp = ft_math(cfg, timelock1);

cfg.value = pi;

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1);

cfg.operation = 'subtract';
tmp = ft_math(cfg, timelock1);

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1);

cfg.operation = 'divide';
tmp = ft_math(cfg, timelock1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with two inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1, timelock2);

cfg.operation = 'subtract';
tmp = ft_math(cfg, timelock1, timelock2);

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1, timelock2);

cfg.operation = 'divide';
tmp = ft_math(cfg, timelock1, timelock2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with more than two inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1, timelock2, timelock1, timelock2);

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1, timelock2, timelock1, timelock2);



