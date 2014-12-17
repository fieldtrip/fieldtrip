function test_bug2785

% WALLTIME 00:10:00
% MEM 1gb

%% create some test data
data = [];
data.label = {'Cz'};
data.time = {
  [0]
  [0 1]
  [0 1 2]
  [0 1 2 3]
  [0 1 2 3 4]
  [0 1 2 3 4 5]
  [0 1 2 3 4 5 6]
  };
data.trial =  {
  [1]
  [2 2]
  [3 3 3]
  [4 4 4 4]
  [5 5 5 5 5]
  [6 6 6 6 6 6]
  [7 7 7 7 7 7 7]
  };

cfg = [];
cfg.vartrllength = 2;
avg = ft_timelockanalysis(cfg, data);

assert(avg.dof(1)==7)
assert(avg.dof(2)==6)
assert(avg.dof(3)==5)
assert(avg.dof(4)==4)
assert(avg.dof(5)==3)
assert(avg.dof(6)==2)
assert(avg.dof(7)==1)

cfg = [];
cfg.method = 'within';
grandavg = ft_timelockgrandaverage(cfg, avg, avg, avg);

assert(grandavg.dof(1)==7*3) % this is at t=0
assert(grandavg.dof(2)==6*3) % this is at t=1
assert(grandavg.dof(3)==5*3) % etc.
assert(grandavg.dof(4)==4*3)
assert(grandavg.dof(5)==3*3)
assert(grandavg.dof(6)==2*3)
assert(grandavg.dof(7)==1*3)

cfg = [];
cfg.method = 'within';
cfg.latency = [0 1];
grandavg = ft_timelockgrandaverage(cfg, avg, avg, avg);

assert(grandavg.dof(1)==7*3) % this is at t=0
assert(grandavg.dof(2)==6*3) % this is at t=1

%% the following section is where it fails

cfg = [];
cfg.method = 'within';
cfg.latency = [1 2];
grandavg = ft_timelockgrandaverage(cfg, avg, avg, avg);

assert(grandavg.dof(1)==6*3) % this is now at t=1
assert(grandavg.dof(2)==5*3) % this is now at t=2

%% pull 52 was not complete and was suspected of having problems with different channels

avg1 = [];
avg1.time = 1:10;
avg1.label = {'a', 'b', 'c'};
avg1.dimord = 'chan_time';
avg1.avg = randn(3, 10);
avg1.dof = [
  [1:10]*1
  [1:10]*10
  [1:10]*100
  ];

avg2 = [];
avg2.time = 6:15;
avg2.label = {'b', 'c', 'd'};
avg2.dimord = 'chan_time';
avg2.avg = randn(3, 10);
avg2.dof = [
  [6:15]*10
  [6:15]*100
  [6:15]*1000
  ];

cfg = [];
cfg.method = 'within';
grandavg = ft_timelockgrandaverage(cfg, avg1, avg2);

assert(grandavg.dof(  1)== 120);
assert(grandavg.dof(end)==2000);
