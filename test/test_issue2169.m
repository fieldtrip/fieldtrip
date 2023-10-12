function test_issue2169

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_apply_montage
% DATA no

%%

data0 = {};
data0.label = {
  '1'
  '2'
  '3'
  '4'
  '5'
  '6'
  '7'
  '8'
  '9'
  };
data0.trial{1} = [
  1 1 1 1 1 1 1 1 1 1
  2 2 2 2 2 2 2 2 2 2
  3 3 3 3 3 3 3 3 3 3
  4 4 4 4 4 4 4 4 4 4
  5 5 5 5 5 5 5 5 5 5
  6 6 6 6 6 6 6 6 6 6
  7 7 7 7 7 7 7 7 7 7
  8 8 8 8 8 8 8 8 8 8
  9 9 9 9 9 9 9 9 9 9 
  ];
data0.time{1} = 1:10;

data1 = data0;
data1.trial{1}(1,10) = nan;

data2 = data0;
data2.trial{1}(2,10) = nan;

data3 = data0;
data3.trial{1}(3,10) = nan;

%%

cfg = [];
cfg.montage.labelold = {'1', '2', '3'};
cfg.montage.labelnew = {'1', '2', '1-2'};
cfg.montage.tra = [
  1  0 0
  0  1 0
  1 -1 0
  ];
reref0 = ft_preprocessing(cfg, data0);
reref1 = ft_preprocessing(cfg, data1);
reref2 = ft_preprocessing(cfg, data2);
reref3 = ft_preprocessing(cfg, data3);

% check that the nans don't spread

assert(~isnan(reref0.trial{1}(1,10)))
assert(~isnan(reref0.trial{1}(2,10)))
assert(~isnan(reref0.trial{1}(3,10)))

assert( isnan(reref1.trial{1}(1,10)))
assert(~isnan(reref1.trial{1}(2,10)))
assert( isnan(reref1.trial{1}(3,10)))

assert(~isnan(reref2.trial{1}(1,10)))
assert( isnan(reref2.trial{1}(2,10)))
assert( isnan(reref2.trial{1}(3,10)))

assert(~isnan(reref3.trial{1}(1,10)))
assert(~isnan(reref3.trial{1}(2,10)))
assert(~isnan(reref3.trial{1}(3,10)))

%%
% repeat the same thing as above, but now for single precision data

data0.trial{1} = single(data0.trial{1});
data1.trial{1} = single(data1.trial{1});
data2.trial{1} = single(data2.trial{1});
data3.trial{1} = single(data3.trial{1});

cfg = [];
cfg.montage.labelold = {'1', '2', '3'};
cfg.montage.labelnew = {'1', '2', '1-2'};
cfg.montage.tra = [
  1  0 0
  0  1 0
  1 -1 0
  ];
reref0 = ft_preprocessing(cfg, data0);
reref1 = ft_preprocessing(cfg, data1);
reref2 = ft_preprocessing(cfg, data2);
reref3 = ft_preprocessing(cfg, data3);

% check that the nans don't spread

assert(~isnan(reref0.trial{1}(1,10)))
assert(~isnan(reref0.trial{1}(2,10)))
assert(~isnan(reref0.trial{1}(3,10)))

assert( isnan(reref1.trial{1}(1,10)))
assert(~isnan(reref1.trial{1}(2,10)))
assert( isnan(reref1.trial{1}(3,10)))

assert(~isnan(reref2.trial{1}(1,10)))
assert( isnan(reref2.trial{1}(2,10)))
assert( isnan(reref2.trial{1}(3,10)))

assert(~isnan(reref3.trial{1}(1,10)))
assert(~isnan(reref3.trial{1}(2,10)))
assert(~isnan(reref3.trial{1}(3,10)))

