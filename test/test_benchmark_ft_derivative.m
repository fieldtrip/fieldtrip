function test_benchmark_ft_derivative

% MEM 1gb
% WALLTIME 00:10:00

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these are the data specific parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_array = [0 1 8 64 128];       % number of channels
n_array = [0 10 100 500 1000];  % number of samples
niter   = 10;                   % number of iterations with the same parameter/variable set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these are the function specific parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

funname = 'ft_preproc_derivative';

clear argname 
argname{1} = 'order';

clear argval 
argval{1} = nan; % see below

% use varying order
argval{1} = 1;
benchmark(funname, argname, argval, m_array, n_array, niter, 'feedback', 'table', 'tableheader', true, 'tabledata', true)
argval{1} = 2;
benchmark(funname, argname, argval, m_array, n_array, niter, 'feedback', 'table', 'tableheader', false, 'tabledata', true)
argval{1} = 3;
benchmark(funname, argname, argval, m_array, n_array, niter, 'feedback', 'table', 'tableheader', false, 'tabledata', true)
