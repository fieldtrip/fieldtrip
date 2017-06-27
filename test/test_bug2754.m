function test_bug2754

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_selectdata

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2754.mat'));

tbeg = 0.0050;
tend = 0.0200;

% fix the data
testRAW = testERF;
testERF = ft_timelockanalysis([], testRAW);

%%

cfg = [];
cfg.avgovertime = 'yes';
avgERF = ft_selectdata(cfg, testERF);

assert(isnumeric(avgERF.time));
assert(length(avgERF.time)==1);

%%

cfg = [];
cfg.latency = [tbeg tend];
avgERF = ft_selectdata(cfg, testERF);

assert(isnumeric(avgERF.time));
assert(length(avgERF.time)==31);
assert(size(avgERF.avg,1)==359 && size(avgERF.avg,2)==31);


%%

cfg = [];
cfg.avgovertime = 'yes';
avgRAW = ft_selectdata(cfg, testRAW);

assert(iscell(avgRAW.time));
assert(size(avgRAW.trial{1},1)==359 && size(avgRAW.trial{1},2)==1 );
assert(size(avgRAW.time{1} ,1)==1   && size(avgRAW.time{1} ,2)==1 );

%%

cfg = [];
cfg.latency = [tbeg tend];
selRAW = ft_selectdata(cfg, testRAW);

assert(iscell(selRAW.time));
assert(length(selRAW.time)==2 && length(selRAW.trial)==2);
assert(size(selRAW.trial{1},1)==359 && size(selRAW.trial{1},2)==31 );
assert(size(selRAW.time{1} ,1)==1   && size(selRAW.time{1} ,2)==31 );

%%

cfg = [];
cfg.latency = [tbeg tend];
cfg.avgovertime = 'yes';
selRAW = ft_selectdata(cfg, testRAW);

assert(iscell(selRAW.time));
assert(length(selRAW.time)==2 && length(selRAW.trial)==2);
assert(size(selRAW.trial{1},1)==359 && size(selRAW.trial{1},2)==1 );
assert(size(selRAW.time{1} ,1)==1   && size(selRAW.time{1} ,2)==1  );

%%

cfg = [];
cfg.latency = [tbeg tend];
cfg.avgovertime = 'yes';
cfg.trials = 1;
selRAW = ft_selectdata(cfg, testRAW);

assert(iscell(selRAW.time));
assert(length(selRAW.time)==1 && length(selRAW.trial)==1);
assert(size(selRAW.trial{1},1)==359 && size(selRAW.trial{1},2)==1 );
assert(size(selRAW.time{1} ,1)==1   && size(selRAW.time{1} ,2)==1  );

cfg.trials = 1;
selRAW1 = ft_selectdata(cfg, testRAW);
cfg.trials = 2;
selRAW2 = ft_selectdata(cfg, testRAW);
assert(selRAW1.trial{1}(1)~=selRAW2.trial{1}(1)); % values should be different

