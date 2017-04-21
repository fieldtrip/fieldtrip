function test_bug3287

% WALLTIME 00:10:00
% MEM 1gb

%%

timelock1 = [];
timelock1.label = {'1'};
timelock1.time  = 1:10;
timelock1.avg   = 1*ones(1,10);
timelock1.trial = 1*ones(3,1,10);

timelock2 = [];
timelock2.label = {'2'};
timelock2.time  = 1:10;
timelock2.avg   = 2*ones(1,10);
timelock2.trial = 2*ones(3,1,10);

timelock3 = [];
timelock3.label = {'1'};
timelock3.time  = 11:20;
timelock3.avg   = 3*ones(1,10);
timelock3.trial = 3*ones(3,1,10);

%%

cfg = [];
append11 = ft_appendtimelock(cfg, timelock1, timelock1); % appenddim = rpt
assert(isequal(size(append11.avg), [2, 1, 10]));
assert(isequal(size(append11.trial), [6, 1, 10]));

cfg = [];
append12 = ft_appendtimelock(cfg, timelock1, timelock2); % appenddim = chan
assert(isequal(size(append12.avg), [2, 10]));
assert(isequal(size(append12.trial), [3, 2, 10]));

try
  cfg = [];
  append12 = ft_appendtimelock(cfg, rmfield(timelock1, 'avg'), rmfield(timelock2, 'trial'));
  catchflag = false;
catch
  catchflag = true;
end
assert(catchflag);

cfg = [];
append13 = ft_appendtimelock(cfg, timelock1, timelock3); % appenddim = time


