function test_bug2597

% WALLTIME 00:10:00
% MEM 1500mb

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
assert(isequal(datanew.label,{'b';'a';'c';'d'}));
assert(isequal(datanew2.label,{'b';'a';'c';'d'}));
assert(isequaln(datanew.avg(:,1),[2 1 3 nan]'));
assert(isequaln(datanew2.avg(:,1),[2 nan 3 4]'));

data2        = [];
data2.label  = {'b';'c';'e'};
data2.avg    = repmat([2 3 5]',[1 2]);
data2.dimord = 'chan_time';
data2.time   = [1 2];

data3        = [];
data3.label  = {'d';'a';'c'};
data3.avg    = repmat([4 1 3]',[1 2]);
data3.dimord = 'chan_time';
data3.time   = [1 2];

cfg.select = 'intersect';
[datanew, datanew2, datanew3] = ft_selectdata(cfg, data, data2, data3);
assert(isequal(datanew.label, {'c'}));
assert(isequal(datanew2.label, {'c'}));
assert(isequal(datanew3.label, {'c'}));
assert(isequaln(datanew.avg(:,1), 3));
assert(isequaln(datanew2.avg(:,1), 3));
assert(isequaln(datanew3.avg(:,1), 3));

cfg.select = 'union';
[datanew, datanew2, datanew3] = ft_selectdata(cfg, data, data2, data3);
assert(isequal(datanew.label, {'b';'a';'c';'e';'d'}));
assert(isequal(datanew2.label, {'b';'a';'c';'e';'d'}));
assert(isequal(datanew3.label, {'b';'a';'c';'e';'d'}));
assert(isequaln(datanew.avg(:,1), [2 1 3 nan nan]'));
assert(isequaln(datanew2.avg(:,1), [2 nan 3 5 nan]'));
assert(isequaln(datanew3.avg(:,1), [nan 1 3 nan 4]'));
