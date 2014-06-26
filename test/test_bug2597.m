function test_bug2597

% WALLTIME 00:10:00
% MEM 150mb

% TEST test_bug2597
% TEST ft_selectdata

% this function tests the functionality of ft_selectdata to preserve the
% ordering of the channels according to the first input argument (thus
% overruling alphabetization, as well as any order in the cfg).

% make some data
data        = [];
data.label  = {'b';'a';'c'};
data.avg    = repmat([2 1 3]',[1 2]);
data.dimord = 'chan_time';
data.time   = [1 2];

datanew = ft_selectdata([], data);
assert(isequal(datanew.label,{'b';'a';'c'}));
assert(isequal(datanew.avg(:,1),[2 1 3]'));

data2        = [];
data2.label  = {'a';'b';'c'};
data2.avg    = repmat([1 2 3]',[1 2]);
data2.dimord = 'chan_time';
data2.time   = [1 2];

[datanew, datanew2] = ft_selectdata([], data, data2);
assert(isequal(datanew.label,{'b';'a';'c'}));
assert(isequal(datanew2.label,datanew.label));

data2        = [];
data2.label  = {'b';'c';'d'};
data2.avg    = repmat([2 3 4]',[1 2]);
data2.dimord = 'chan_time';
data2.time   = [1 2];
[datanew, datanew2] = ft_selectdata([], data, data2);
assert(isequal(datanew.label,{'b';'c'}));
assert(isequal(datanew2.label,{'b';'c'}));
assert(isequal(datanew.avg(:,1),[2 3]'));
assert(isequal(datanew2.avg(:,1),[2 3]'));

cfg.select = 'union';
[datanew, datanew2] = ft_selectdata(cfg, data, data2);
