function test_ft_selectdata_union

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_selectdata

data1        = [];
data1.label  = {'1' '2'};
data1.time   = 1:5;
data1.dimord = 'chan_time';
data1.avg    = 1*ones(2,5);

data2        = [];
data2.label  = {'1' '2'};
data2.time   = 2:6;         % NOTE THE DIFFERENCE
data2.dimord = 'chan_time';
data2.avg    = 2*ones(2,5);

data3        = [];
data3.label  = {'2' '3'};   % NOTE THE DIFFERENCE
data3.time   = 1:5;
data3.dimord = 'chan_time';
data3.avg    = 3*ones(2,5);

data4        = [];
data4.label  = {'2' '3'};   % NOTE THE DIFFERENCE
data4.time   = 2:6;         % NOTE THE DIFFERENCE
data4.dimord = 'chan_time';
data4.avg    = 4*ones(2,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first with intersect

cfg = [];
cfg.channel = 'all';
cfg.latency = 'all';
cfg.select = 'intersect';

[a, b] = ft_selectdata(cfg, data1, data2);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==2);
assert(length(a.time)==4);

[a, b] = ft_selectdata(cfg, data1, data3);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==1);
assert(length(a.time)==5);

[a, b] = ft_selectdata(cfg, data1, data4);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==1);
assert(length(a.time)==4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% then with union

cfg = [];
cfg.channel = 'all';
cfg.latency = 'all';
cfg.select = 'union';

[a, b] = ft_selectdata(cfg, data1, data2);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==2);
assert(length(a.time)==6);

[a, b] = ft_selectdata(cfg, data1, data3);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==3);
assert(length(a.time)==5);

[a, b] = ft_selectdata(cfg, data1, data4);

assert(isequal(a.label, b.label));
assert(isequal(a.time, b.time));
assert(length(a.label)==3);
assert(length(a.time)==6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now with three dimensions

clear all

data1        = [];
data1.label  = {'1' '2'};
data1.time   = 1:5;
data1.freq   = 1:3;
data1.dimord = 'chan_freq_time';
data1.powspctrm = 1*ones(2,3,5);

data2        = [];
data2.label  = {'2' '3'};
data2.time   = 2:6;
data2.freq   = 2:4;
data2.dimord = 'chan_freq_time';
data2.powspctrm = 2*ones(2,3,5);

cfg = [];
cfg.select = 'union';
[a, b] = ft_selectdata(cfg, data1, data2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now with two source inputs

if false
  % this does not yet work on 11 April 2014
  data1 = [];
  data1.pos = randn(5,3);
  data1.pow = 1*ones(5,1);
  
  data2 = [];
  data2.pos = randn(5,3);
  data2.pow = 2*ones(5,1);
  
  cfg = [];
  cfg.select = 'union';
  [a, b] = ft_selectdata(cfg, data1, data2);
  
end





