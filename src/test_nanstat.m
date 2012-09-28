function test_suite = test_nanstat
% TEST nansum nanmean nanvar nanstd
%
% Test conformance of FieldTrip's nansum, nanmean, nanvar & nanstd
% functions with MATLABs versions in the statistics toolbox.

initTestSuite;  % for xUnit

function setup
fprintf('Extending path...\n');
addpath ../src

function teardown
fprintf('Cleaning up path...\n');
rmpath ../src


function test_datatypes
% test different data types
X = {};
X{end+1} = [true, true, false];
X{end+1} = 'a string';
X{end+1} = uint8(1:25);
%X{end+1} = int64(1:25);    % matlab 2010a can't sum 64bit ints :/
X{end+1} = complex(rand(10), rand(10));  % test complex numbers
X{end+1} = cast(X{end}, 'single');  % test single precision complex numbers
X{end+1} = 0;
X{end+1} = inf;
X{end+1} = -inf;

for i = 1:length(X)
  x = X{i};
  class(x), class(nansum(x)), class(sum(x))
  fprintf('norm = %.3g', norm(nansum(x) - sum(x)));
  assertEqual(nansum(x), sum(x));
  %assertEqual(nanmean(x), mean(x));
  
%   if (length(strmatch(class(x), {'uint8', 'int64'})) > 0)
%     fprintf('Skipping type %s for nanvar & nanstd\n', class(x));
%   else
%     if isinf(x)
%       % Special case since at least MATLAB Version 7.11.0.584 (R2010b)
%       % seems give different results regular variants of the nan
%       % statistics:
%       assertEqual(nanvar(x), 0);
%       assertEqual(nanstd(x), 0);
%     else
%       %fprintf('Testing type %s for nanvar & nanstd\n', class(x));
%       assertEqual(nanvar(x), var(x));
%       assertEqual(nanstd(x), std(x));  
%     end
%   end
end


function test_variance_precision
% test for inaccuracy caused by running estimates
signal = rand(1000, 1);

for offset = logspace(8, 20, 4)
  x = signal + offset;
  sig2_true = var(x);  % note that var also suffers from numerical imprecision.
  sig2 = nanvar(x);
  fprintf('Adding offset %.2g: var()=%.4f, nanvar()=%.4f.\n', ...
    offset, sig2_true, sig2);
  assert(abs(sig2 - sig2_true) < eps, sprintf(...
    ['Numerical imprecision detected in nanvar: '...
    '%.4g != %.4g at offset %.2g.'], sig2, sig2_true, offset));
end


function test_nansum
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN


assert(nansum([1, 2, NaN, 3]) == 6);  % vector -> scalar
assertElementsAlmostEqual(nansum(X), [7, 6, 0]);
assertElementsAlmostEqual(nansum(reshape(X, [1 3 3])), ...
  reshape([7, 6, 0], [1 1 3]));
assertElementsAlmostEqual(nansum(X, 1), [7, 6, 0]);
assertElementsAlmostEqual(nansum(X, 2), [1, 8, 4]');

Y = nansum(X, 3); % same size, but nans -> 0.
assertElementsAlmostEqual(Y(isnan(X)), [0, 0, 0, 0, 0]');
assertElementsAlmostEqual(Y(~isnan(X)), X(~isnan(X)));


function test_nanmean
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


function test_nanvar
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


function test_nanstd
X = magic(3);
X([1 6:9]) = repmat(NaN,1,5);
% X =
%    NaN     1   NaN
%      3     5   NaN
%      4   NaN   NaN

assertElementsAlmostEqual(nanstd(X), sqrt(nanvar(X)));
assertElementsAlmostEqual(nanstd(X, 1), sqrt(nanvar(X, 1)));
assertElementsAlmostEqual(nanstd(X * NaN, 1), sqrt(nanvar(X * NaN, 1)));
