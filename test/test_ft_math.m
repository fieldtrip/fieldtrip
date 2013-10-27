function test_ft_math

% MEM 1gb
% WALLTIME 00:03:04

% TEST test_ft_math
% TEST ft_math

% create some data
timelock1.label  = {'chan1'; 'chan2'};
timelock1.time   = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg    = ones(2,5);
timelock1.cfg    = struct([]);

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
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.value = pi;

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'subtract';
tmp = ft_math(cfg, timelock1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'divide';
tmp = ft_math(cfg, timelock1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with two inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'subtract';
tmp = ft_math(cfg, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'divide';
tmp = ft_math(cfg, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with more than two inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'add';
tmp = ft_math(cfg, timelock1, timelock2, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1, timelock2, timelock1, timelock2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

