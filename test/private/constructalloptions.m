function optarg = constructalloptions(arg)

% CONSTRUCTALLOPTIONS allows you to create an extensive set of key-value
% pairs that can be used for testing a function.
%
% Use as
%   optarg = test_constructoptarg(arg)
% where
%   arg(1).name  = string, name of the first option
%   arg(1).value = cell-array with all possible values for the first option
%   arg(2).name  = string, name of the second option
%   arg(2).value = cell-array with all possible values for the second option
%   ...
%
% See also CFG2KEYVAL, KEYVAL2CFG

optarg = {};
numopt = numel(arg);
numpre = 0;

for i=numopt:-1:1
  numval = numel(arg(i).value);
  fprintf('adding option %d with %d values, ', i, numval);
  if isempty(optarg)
    numpre = 1;
    optarg = cell(numval, numopt*2);
  else
    numpre = size(optarg,1);
    optarg = repmat(optarg,numval,1);
  end
  for j=1:numval
    blocksel = numpre*(j-1) + (1:numpre);
    optarg(blocksel,2*i-1) = {arg(i).name};
    optarg(blocksel,2*i  ) = {arg(i).value{j}};
  end % for numvals
  fprintf('%d in total\n', size(optarg,1));
end % for numopt


% for i=1:numopt
%   numval = numel(arg(i).value);
%   fprintf('adding option %d with %d values to the previous %d configurations\n', i, numval, size(optarg,1));
%   for j=1:max(size(optarg,1),1)
%     if isempty(optarg)
%       thisopt = cell(numval, numopt*2);
%     else
%       thisopt = repmat(optarg(j,:),numval,1);
%     end
%     thisopt(:, 2*i-1) = repmat({arg(i).name},numval,1);
%     thisopt(:, 2*i  ) = arg(i).value(:);
%     optarg = cat(1, optarg, thisopt);
%   end
% end

