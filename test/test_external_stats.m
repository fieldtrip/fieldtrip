function test_external_stats

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY hanning

filelist = {'betacdf'
            'betainv'
            'betapdf'
            'binocdf'
            'binopdf'
            'finv'
            'mvnrnd'
            'nanmax'
            'nanmean'
            'nanmin'
            'nanstd'
            'nansum'
            'nanvar'
            'range'
            'tcdf'
            'tinv'};

[ftver, ftpath] = ft_version;
restoredefaultpath

% ensure that the the 'compat' (i.e. external/signal) is used
global ft_default
ft_default.toolbox.signal = 'compat';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  assert(exist(filelist{k}, 'file')==2||exist(filelist{k}, 'file')==3);
  funhandle = str2func(filelist{k});

  fprintf('testing the functionality of %s\n', filelist{k});
  switch filelist{k}
    case 'betacdf'
      x = [-1 0 0.5 1 2];
      y = [0 0 0.75 1 1];
      assert(isequal(betacdf(x, ones(1,5), 2*ones(1,5)), y));
      assert(isequal(betacdf(x, 1, 2*ones(1,5)), y));
      assert(isequal(betacdf(x, ones(1,5), 2), y));
      assert(isequaln(betacdf(x, [0 1 NaN 1 1], 2), [NaN 0 NaN 1 1]));
      assert(isequaln(betacdf(x, 1, 2*[0 1 NaN 1 1]), [NaN 0 NaN 1 1]));
      assert(isequaln(betacdf([x(1:2) NaN x(4:5)], 1, 2), [y(1:2) NaN y(4:5)]));

      % Test class of input preserved
      assert(isequaln(betacdf([x, NaN], 1, 2), [y, NaN]));
      assert(isequaln(betacdf(single([x, NaN]), 1, 2), single([y, NaN])));
      assert(isequaln(betacdf([x, NaN], single(1), 2), single([y, NaN])));
      assert(isequaln(betacdf([x, NaN], 1, single(2)), single([y, NaN])));
      
      input.betacdf   = {rand(10,20)*10 1 2};
    case 'betainv'
      input.betainv   = {rand(10,20) 1 2};
    case 'betapdf'
      x = [-1 0 0.5 1 2];
      y = [0 2 1 0 0];
      assert(isequal(betapdf(x, ones(1,5), 2*ones(1,5)), y));
      assert(isequal(betapdf(x, 1, 2*ones(1,5)), y));
      assert(isequal(betapdf(x, ones(1,5), 2), y));
      assert(isequaln(betapdf(x, [0 NaN 1 1 1], 2), [NaN NaN y(3:5)]));
      assert(isequaln(betapdf(x, 1, 2*[0 NaN 1 1 1]), [NaN NaN y(3:5)]));
      assert(isequaln(betapdf([x, NaN], 1, 2), [y, NaN]));

      % Test class of input preserved
      assert(isequaln(betapdf (single([x, NaN]), 1, 2), single([y, NaN])));
      assert(isequaln(betapdf ([x, NaN], single(1), 2), single([y, NaN])));
      assert(isequaln(betapdf ([x, NaN], 1, single(2)), single([y, NaN])));

      % Beta (1/2,1/2) == arcsine distribution
      x = rand (10,1);
      y = 1./(pi * sqrt (x.*(1-x)));
      assert(isalmostequal(betapdf(x, 1/2, 1/2), y, 'abstol', 50*eps));

      % Test large input values to betapdf
      assert(isalmostequal(betapdf(0.5, 1000, 1000), 35.678, 'abstol', 1e-3));

      input.betapdf   = {rand(10,20)*10 1 2};
    case 'binocdf'
      input.(filelist{k})   = {4 6 rand(10,20)};
    case 'binopdf'
      input.(filelist{k})   = {4 6 rand(10,20)};
    case 'finv'
      input.(filelist{k})   = {randn(10,20) 30 2};
    case 'mvnrnd'
      input.(filelist{k})   = {randn(10,1) abs(randn(10,1))};
    case 'nanmax'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7) = nan;
      input.(filelist{k})   = {tmp [] 1};
    case 'nanmean'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7)  = nan;
      input.(filelist{k})   = {tmp 1};
    case 'nanmin'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7)  = nan;
      input.(filelist{k})   = {tmp [] 2};
    case 'nanstd'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7)  = nan;
      input.(filelist{k})   = {tmp [] 2};
    case 'nansum'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7)  = nan;
      input.(filelist{k})   = {tmp 1};
    case 'nanvar'
      tmp = randn(10,20);
      tmp(rand(10,20)>0.7)  = nan;
      input.(filelist{k})   = {tmp [] 1};
    case 'range'
      input.(filelist{k})   = {randn(10,20) 1};
    case 'tcdf'
      input.(filelist{k})   = {randn(10,20) 75};
    case 'tinv'
      input.(filelist{k})   = {rand(10,20) 75};
    otherwise
      ft_error('function %s is not part of the official external/stats directory', filelist{k});
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the external/stats output with the matlab version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:numel(filelist)
  try
    % this only works for the functions that create a window
    funhandle = str2func(filelist{k});
    output1.(filelist{k}) = funhandle(input.(filelist{k}){:});
  end
end

restoredefaultpath

% ensure that the the matlab version is used
ft_default.toolbox.signal = 'matlab';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  try
    % this only works for the functions that create a window
    funhandle = str2func(filelist{k});
    output2.(filelist{k}) = funhandle(input.(filelist{k}){:});
  end
end

fn1 = fieldnames(output1);
fn2 = fieldnames(output2);
assert(isequal(sort(fn1), sort(fn2)));
for k = 1:numel(fn1)
  assert(isalmostequal(output1.(fn1{k}), output2.(fn1{k}), 'abstol', 10*eps));
end
