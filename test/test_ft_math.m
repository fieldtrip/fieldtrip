function test_ft_math

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_math

% create some test data
raw1 = [];
raw1.time = {[1:10], [1:10], [1:10]};
raw1.trial = {ones(2,10), ones(2,10), ones(2,10)};
raw1.label = {'chan01';'chan02'};
raw1.trialinfo = rand(3,4);

timelock1.label  = {'chan1'; 'chan2'};
timelock1.time   = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg    = ones(2,5);
timelock1.cfg    = struct([]);

timelock2 = timelock1;
timelock2.avg  = ones(2,5)*2;

timelock3 = timelock1;
timelock3.var = ones(2,5);

source1 = [];
source1.pos = randn(10,3);
source1.pow = randn(10,1);
source1.powdimord = 'pos';
source1.mom = cell(10,1);
for i=1:10
  source1.mom{i} = ones(3, 20)*1;
end
source1.momdimord = '{pos}_ori_time';
source1.time = 1:20;

source2 = source1;
for i=1:10
  source2.mom{i} = ones(3, 20)*2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with one raw input

cfg=[];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter = 'trial';

cfg.operation = 'log10';
tmp = ft_math(cfg, raw1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');

cfg.scalar = pi;

cfg.operation = 'add';
tmp = ft_math(cfg, raw1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');

cfg.operation = 'subtract';
tmp = ft_math(cfg, raw1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, raw1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');

cfg.operation = 'divide';
tmp = ft_math(cfg, raw1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with one timelock input

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'avg';

cfg.operation = 'log10';
tmp = ft_math(cfg, timelock1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');

cfg.scalar = pi;

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
% do operation with one timelock input, multiple parameters

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = {'avg', 'var'};

cfg.operation = 'log10';
tmp = ft_math(cfg, timelock3);
assert(isfield(tmp, 'avg') && isfield(tmp, 'var'), 'the output parameter is missing');
assert(isfield(tmp, 'dimord'), 'the output dimord is missing');
assert(isequal(tmp.var, log10(timelock3.var)));
assert(isequal(tmp.avg, log10(timelock3.avg)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with one source input

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'mom'; % note that this is {pos}_ori_time

cfg.operation = 'log10';
tmp = ft_math(cfg, source1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.scalar = pi;

cfg.operation = 'add';
tmp = ft_math(cfg, source1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'subtract';
tmp = ft_math(cfg, source1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, source1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'divide';
tmp = ft_math(cfg, source1);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with two timelock inputs

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
% do operation with two source inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'mom';

cfg.operation = 'add';
tmp = ft_math(cfg, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'subtract';
tmp = ft_math(cfg, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'divide';
tmp = ft_math(cfg, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with more than two timelock inputs

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do operation with more than two source inputs

cfg = [];
cfg.showcallinfo = 'no';
cfg.trackconfig  = 'no';
cfg.parameter    = 'mom';

cfg.operation = 'add';
tmp = ft_math(cfg, source1, source2, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');

cfg.operation = 'multiply';
tmp = ft_math(cfg, source1, source2, source1, source2);
assert(isfield(tmp, cfg.parameter), 'the output parameter is missing');
assert(isfield(tmp, [cfg.parameter 'dimord']), 'the output dimord is missing');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the numerical output of the operation

cfg = [];
cfg.parameter = 'trial';

cfg.operation = 'log10';
tmp = ft_math(cfg, raw1);
assert(tmp.trial{1}(1)==0);

cfg.operation = 'multiply';
cfg.scalar = -1;
tmp = ft_math(cfg, raw1);
assert(tmp.trial{1}(1)==-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idem for a timelock input

cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'add';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==3);

cfg.operation = 'subtract';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==-1);

cfg.operation = 'divide';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==1/2);

cfg.operation = 'multiply';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==2);

cfg.operation = 'log10';
tmp = ft_math(cfg, timelock1);
assert(tmp.avg(1)==0);

cfg.operation = 'log(x1)';
tmp = ft_math(cfg, timelock1);
assert(tmp.avg(1)==0);

cfg.operation = '(x1-x2)/(x1+x2)';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==-1/3);

cfg.scalar = 2;
cfg.operation = '(x1+x2)^s';
tmp = ft_math(cfg, timelock1, timelock2);
assert(tmp.avg(1)==9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idem for a source structure with a cell array

cfg = [];
cfg.parameter = 'mom';
cfg.operation = 'add';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==3);

cfg.operation = 'subtract';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==-1);

cfg.operation = 'divide';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==1/2);

cfg.operation = 'multiply';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==2);

cfg.operation = 'log10';
tmp = ft_math(cfg, source1);
assert(tmp.mom{1}(1)==0);

cfg.operation = 'log(x1)';
tmp = ft_math(cfg, source1);
assert(tmp.mom{1}(1)==0);

cfg.operation = '(x1-x2)/(x1+x2)';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==-1/3);

cfg.scalar = 2;
cfg.operation = '(x1+x2)^s';
tmp = ft_math(cfg, source1, source2);
assert(tmp.mom{1}(1)==9);