function test_nanstat

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_nanstat nansum nanmean nanstd nanvar nanvar_base

% Test the conformance of FieldTrip's nansum, nanmean, nanvar and nanstd
% functions with MATLABs versions in the statistics toolbox.

% this includes assertVectorsAlmostEqual etcetera
ft_hastoolbox('xunit', 1);

[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'src'));

% test different data types
X = {};
X{end+1} = [true, true, false];
X{end+1} = 'a string';
X{end+1} = single(1:25);
X{end+1} = double(1:25);
X{end+1} = int8(1:25);
X{end+1} = uint8(1:25);
X{end+1} = int16(1:25);
X{end+1} = uint16(1:25);
X{end+1} = int32(1:25);
X{end+1} = uint32(1:25);
%X{end+1} = int64(1:25);                  % MATLAB 2010a and older can't sum 64bit ints :/
%X{end+1} = uint64(1:25);                 % MATLAB 2010a and older can't sum 64bit ints :/
X{end+1} = complex(rand(10), rand(10));   % test double precision complex numbers
X{end+1} = cast(X{end}, 'single');        % test single precision complex numbers
X{end+1} = 0; 
X{end+1} = inf;
X{end+1} = -inf;
X{end+1} = [1 2 inf]; % FIXME: what is expected here?

for i = 1:length(X)
  x = X{i};
  assertEqual(nansum(x), sum(x));
  assertEqual(nanmean(x), mean(x));
  
  if isinteger(x)
    fprintf('Skipping type %s for nanvar & nanstd\n', class(x));
    % MATLAB Version 2011b does not support var and std for integer input
    % types.
  elseif any(isinf(x))
    fprintf('Skipping nanvar & nanstd for inf\n');
    % Special case, since at least MATLAB Version 7.11.0.584 (R2010b) seems
    % give different results regular variants of the nan statistics.
    % Results seem rather arbitrary, so we don't emulate. Skipping.
  elseif x==0
    fprintf('Skipping nanvar & nanstd for 0\n');
    % test case fails, var(0) results in 0, nanvar(0) results in NaN. This specific case could be catched in 
    % our mex implementation but right now (24-09-2013, roevdmei+jansch) we decided to skip it in favor
    % of having the test function run to detect more serious errors
  else
    fprintf('Testing type %s for nanvar & nanstd\n', class(x));
    assertElementsAlmostEqual(nanvar(x), var(x));
    assertElementsAlmostEqual(nanstd(x), std(x));
  end
end

% test for inaccuracy caused by running estimates
signal = rand(1000, 1);

for offset = logspace(8, 12, 4)
  x = signal + offset;
  sig2_true = var(x);  % note that var also suffers from (serious) numerical imprecision.
  sig2 = nanvar(x);
  fprintf('Adding offset %.2g: var()=%.4f, nanvar()=%.4f.\n', offset, sig2_true, sig2);
  assert(isalmostequal(sig2_true, sig2, 'reltol', 1e-3));
end

% test nansum
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN

assert(nansum([1, 2, NaN, 3]) == 6);  % vector -> scalar
assertElementsAlmostEqual(nansum(X), [7, 6, 0]);
assertElementsAlmostEqual(nansum(reshape(X, [1 3 3])), reshape([7, 6, 0], [1 1 3]));
assertElementsAlmostEqual(nansum(X, 1), [7, 6, 0]);
assertElementsAlmostEqual(nansum(X, 2), [1, 8, 4]');

Y = nansum(X, 3); % same size, but nans -> 0.
assertElementsAlmostEqual(Y(isnan(X)), [0, 0, 0, 0, 0]');
assertElementsAlmostEqual(Y(~isnan(X)), X(~isnan(X)));

% test nanmean
% Extended example from nanmean documentation:
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN

assertElementsAlmostEqual(nanmean(X), [3.5, 3, NaN]);
assertElementsAlmostEqual(nanmean(X, 1), [3.5, 3, NaN]);
assertElementsAlmostEqual(nanmean(X, 2), [1, 4, 4]');
assertElementsAlmostEqual(nanmean(X, 3), X);
assertElementsAlmostEqual(nanmean(X * NaN), [NaN, NaN, NaN])

% Test higher order dimensions:
X = randn(2, 3, 5, 7);
for dim = 1:4
  assertElementsAlmostEqual(nanmean(X, dim), mean(X, dim));
end

%{
% the following tests are disabled as our nanvar treats dimensions with all
nans a bit different from mathworks; nanvar([nan nan nan]) = 0 in our case;
nan in mathworks' case. 0 makes more sense imho (and exotic use case
anyway).

Eelke Spaak, 24 oct 2013

% test nanvar
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN

assertElementsAlmostEqual(nanvar(X), [.5, 8, NaN]);
assertElementsAlmostEqual(nanvar(X, 1), [.25, 4, NaN]);
% higher dims throw error
assertElementsAlmostEqual(nanvar(X * NaN), [NaN, NaN, NaN]);

%}

% Test different call signatures & var compatibility
X = randn(2, 3, 5, 7);
for dim = 1:4
  assertElementsAlmostEqual(nanvar(X, 0, dim), var(X, 0, dim));
  assertElementsAlmostEqual(nanvar(X, 1, dim), var(X, 1, dim));
  assertElementsAlmostEqual(nanvar(X, [], dim), var(X, [], dim));
  warning('Also test vector for w!');
end

% test nanstd
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN

assertElementsAlmostEqual(nanstd(X), sqrt(nanvar(X)));
assertElementsAlmostEqual(nanstd(X, 1), sqrt(nanvar(X, 1)));
assertElementsAlmostEqual(nanstd(X * NaN, 1), sqrt(nanvar(X * NaN, 1)));

% Test different call signatures & std compatibility
X = randn(2, 3, 5, 7);
for dim = 1:4
  assertElementsAlmostEqual(nanstd(X, 0, dim),  std(X, 0, dim));
  assertElementsAlmostEqual(nanstd(X, 1, dim),  std(X, 1, dim));
  assertElementsAlmostEqual(nanstd(X, [], dim), std(X, [], dim));
  warning('Also test vector for w!');
end
