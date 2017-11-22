function test_bug3287

% WALLTIME 00:10:00
% MEM 1gb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1 - test the appending of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = [];
data1.label = {'1'};
data1.time{1}  = 1:10;
data1.trial{1} = 1*ones(1,10);
data1.time{2}  = 1:10;
data1.trial{2} = 1*ones(1,10);
data1.time{3}  = 1:10;
data1.trial{3} = 1*ones(1,10);
data1.trialinfo = [
  1 1 
  1 2
  1 3
  ];

data2 = [];
data2.label = {'2'};
data2.time{1}  = 1:10;
data2.trial{1} = 2*ones(1,10);
data2.time{2}  = 1:10;
data2.trial{2} = 2*ones(1,10);
data2.time{3}  = 1:10;
data2.trial{3} = 2*ones(1,10);
data2.trialinfo = [
  2 1 
  2 2
  2 3
  ];

data3 = data1;
data3.time{1}  = 11:20; % shift the time axis relative to data1
data3.time{2}  = 11:20;
data3.time{3}  = 11:20;

%%

cfg = [];
append11 = ft_appenddata(cfg, data1, data1); % appenddim = rpt
assert(isequal(length(append11.trial), 6));
assert(isfield(append11, 'trialinfo'));
assert(isequal(size(append11.trialinfo), [6, 2]));

cfg = [];
append12 = ft_appenddata(cfg, data1, data2); % appenddim = chan
assert(isequal(length(append12.label), 2));

cfg = [];
append13 = ft_appenddata(cfg, data1, data3);
assert(isequal(length(append13.trial), 6));
assert(isfield(append13, 'trialinfo'));
assert(isequal(size(append13.trialinfo), [6, 2]));

try
  cfg = [];
  append23 = ft_appenddata(cfg, data2, data3);
  catchflag = false;
catch
  catchflag = true;
end
assert(catchflag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2 - the same with timelocked data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

cfg = [];
cfg.keeptrials = 'yes';
timelock1 = ft_timelockanalysis(cfg, data1);
timelock2 = ft_timelockanalysis(cfg, data2);
timelock3 = ft_timelockanalysis(cfg, data3);


%%

cfg = [];
append11 = ft_appendtimelock(cfg, timelock1, timelock1); % appenddim = rpt
assert(isequal(size(append11.trial), [6, 1, 10]));
assert(isfield(append11, 'trialinfo'));
assert(isequal(size(append11.trialinfo), [6, 2]));


cfg = [];
append12 = ft_appendtimelock(cfg, timelock1, timelock2); % appenddim = chan
assert(isequal(size(append12.trial), [3, 2, 10]));
assert(~isfield(append12, 'trialinfo'));

try
  cfg = [];
  cfg.parameter = {'avg', 'trial'};
  append12 = ft_appendtimelock(cfg, rmfield(timelock1, 'avg'), rmfield(timelock2, 'trial'));
  catchflag = false;
catch
  catchflag = true;
end
assert(catchflag);

cfg = [];
append13 = ft_appendtimelock(cfg, timelock1, timelock3); % appenddim = time
assert(isfield(append13, 'trialinfo'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3 - shuffle channel labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = [];
data1.label = {'1', '2', '3'};
data1.time{1}  = 1:10;
data1.trial{1} = 1*ones(3,10);
data1.time{2}  = 1:10;
data1.trial{2} = 1*ones(3,10);
data1.time{3}  = 1:10;
data1.trial{3} = 1*ones(3,10);

data2 = [];
data2.label = {'2', '3', '1'};
data2.time{1}  = 1:10;
data2.trial{1} = 1*ones(3,10);
data2.time{2}  = 1:10;
data2.trial{2} = 1*ones(3,10);
data2.time{3}  = 1:10;
data2.trial{3} = 1*ones(3,10);

append12 = ft_appenddata([], data1, data2);
assert(isequal(append12.label, data1.label));

