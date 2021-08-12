function tests = test_ft_connectivity_mutualinformation

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_mutualinformation

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_chan_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan  = 6;
ntime  = 100000;
dat    = randn(nchan, ntime);

%   histmethod = The way that histograms are generated from the data. Possible values
%                are 'eqpop' (default), 'eqspace', 'ceqspace', 'gseqspace'.
%                See the help of the 'binr' function in the ibtb toolbox for more information.
%   numbin     = scalar value. The number of bins used to create the histograms needed for
%                the entropy computations
%   opts       = structure that is passed on to the 'information' function in the ibtb
%                toolbox. See the help of that function for more information.
%   refindx    = scalar value or 'all'. The channel that is used as 'reference channel'.

result = {};
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'eqpop'    , 'numbin', 5);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'eqspace'  , 'numbin', 5);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'ceqspace' , 'numbin', 5);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'gseqspace', 'numbin', 5);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'eqpop'    , 'numbin', 50);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'eqspace'  , 'numbin', 50);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'ceqspace' , 'numbin', 50);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'ibtb', 'histmethod', 'gseqspace', 'numbin', 50);
result{end+1} = ft_connectivity_mutualinformation(dat, 'method', 'gcmi');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
