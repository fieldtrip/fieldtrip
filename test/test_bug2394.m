function test_bug2394

% WALLTIME 00:10:00
% MEM 1gb

% TEST test_bug2394 ft_selectdata ft_selectdata_new ft_selectdata_old

%% create the data

data = [];
data.label = {'a'};
data.time = repmat({1:1000},100,1);
data.trial = repmat({randn(100,1000)},100,1);

%% subselect ft_selectdata_old

%{
% don't check this anymore
datsel = ft_selectdata_old(data, 'rpt', []);
if numel(datsel.trial) == 0;
  oldnum = 'success';
else
  oldnum = 'fail';
end

datsel = ft_selectdata_old(data, 'rpt', false(numel(data.trial),1));
if numel(datsel.trial) == 0;
  oldlog = 'success';
else
  oldlog = 'fail';
end
%}

%% subselect ft_selectdata_new

cfg = [];
cfg.trials = [];
datsel = ft_selectdata_new(cfg, data);
if numel(datsel.trial) == 0;
  newnum = 'success';
else
  newnum = 'fail';
end

cfg = [];
cfg.trials = false(100,1);
datsel = ft_selectdata_new(cfg, data);
if numel(datsel.trial) == 0;
  newlog = 'success';
else
  newlog = 'fail';
end

%% report

%fprintf('ft_selectdata_old, numeric indexing: %s\n', oldnum);
%fprintf('ft_selectdata_old, logical indexing: %s\n', oldlog);
fprintf('ft_selectdata_new, numeric indexing: %s\n', newnum);
fprintf('ft_selectdata_new, logical indexing: %s\n', newlog);

if strfind([newnum newlog], 'fail')
  error('the test failed, see above for details');
end
