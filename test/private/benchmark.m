function benchmark(funname, argname, argval, m_array, n_array, niter, varargin)

% BENCHMARK a given function
%
% Use as
%   benchmark(funname, argname, argval, m_array, n_array, niter, ...)
%
% Optional input arguments should come in key-value pairs and may include
%   feedback    = none, figure, text, table, all
%   tableheader = true, false
%   tabledata   = true, false
%   selection   = 3x2 array with nchans and nsamples to be used for the table

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the optional input arguments
feedback    = keyval('feedback',    varargin); % none, figure, text, table, all
tableheader = keyval('tableheader', varargin); % true, false
tabledata   = keyval('tabledata',   varargin); % true, false
selection   = keyval('selection',   varargin); % 3x2 array with nchans and nsamples to be used for the table

% set the defaults
if isempty(feedback)
  feedback = 'all';
end
if isempty(tableheader)
  tableheader = true;
end
if isempty(tabledata)
  tabledata = true;
end
if isempty(selection)
  selection = [
    8 100
    8 500
    64 500
    ];
end

% convert the function from a string to a handle
funhandle = str2func(funname);

% this will hold the time that all computations took
t_array   = nan(length(m_array), length(n_array));

%  do the actual benchmarking
for m_indx=1:length(m_array)
  for n_indx=1:length(n_array)

    m = m_array(m_indx);
    n = n_array(n_indx);

    if strcmp(feedback, 'table')
      if ~any(selection(:,1)==m & selection(:,2)==n)
        continue
      end
    end

    % create some random data
    dat = randn(m, n);

    elapsed = zeros(1,niter);
    for iteration=1:niter
      tic;
      funhandle(dat, argval{:});
      elapsed(iteration) = toc*1000; % convert from s into ms
    end

    % remember the amount of time spent on the computation for this M and N
    t_array(m_indx, n_indx) = robustmean(elapsed);

    % give some feedback on screen
    if strcmp(feedback, 'text') || strcmp(feedback, 'all')
      fprintf('nchans = %d, nsamples = %d, time = %f ms\n', m, n, t_array(m_indx, n_indx));
    end

  end
end

if strcmp(feedback, 'figure') || strcmp(feedback, 'all')
  %  give some output in a figure
  figure
  surf(n_array, m_array, t_array);
end

if strcmp(feedback, 'table') || strcmp(feedback, 'all')
  %  give some output to screen that can be copied and pasted into the wiki

  m1 = find(m_array==selection(1,1));  % channels
  n1 = find(n_array==selection(1,2));  % samples
  m2 = find(m_array==selection(2,1));  % channels
  n2 = find(n_array==selection(2,2));  % samples
  m3 = find(m_array==selection(3,1));  % channels
  n3 = find(n_array==selection(3,2));  % samples

  if tableheader
    fprintf('^function name and algorithm details ^ %dch x %dsmp ^ %dch x %dsmp ^ %dch x %dsmp ^\n', ...
      m_array(m1), n_array(n1), ...
      m_array(m2), n_array(n2), ...
      m_array(m3), n_array(n3));
  end
  if tabledata
    str = [];
    dum = sprintf('%s;\n', funname);
    str = cat(2, str, dum);
    for i=1:length(argval)
      dum = printstruct(argname{i}, argval{i});
      str = cat(2, str, dum);
    end
    str(str==10) = ' ';
    fprintf('|%s |  %.2f ms  |  %.2f ms  |  %.2f ms  |\n', str, ...
      t_array(m1, n1), ...
      t_array(m2, n2), ...
      t_array(m3, n3));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for robust estimation of mean, removing outliers on both sides
% select the central part of the sorted vector, a quarter of the values is removed from both sides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = robustmean(x)
x    = sort(x);
n    = length(x);
trim = round(0.25*n);
sel  = (trim+1):(n-trim);
y    = mean(x(sel));
