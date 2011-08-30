function output = memtic(action, counter)

% MEMTIC start a MATLAB memory recorder
%
%   MEMTIC and MEMTOC functions work together to measure memory usage.
%   MEMTIC, by itself, saves the current memory footprintt that MEMTOC
%   uses later to measure the memory that was used between the two.
%
% Use as
%   MEMTIC
%   MEMTOC
% to print the estimated memory use on screen, or
%   MEMTIC
%   M = MEMTOC
% to return the estimated memory (in bytes) in variable M, or
%   C = MEMTIC
%   M = MEMTOC(C)
% to specifically estimate the memory use between a well-defined tic/toc pair.
%
% Note that MATLAB uses internal memory allocation, garbage collection, shallow
% copies of variables, and virtual memory. Due to the advanced handling of
% memory for its variables, it is not easy and in certain cases not possible to
% make a reliable and reproducible estimate based on the memory information
% provided by the operating system.
%
% Example: measure the memory increase due to allocating a lot of memory.
% Doing a "clear x" following the allocation and priot to MEMTOC does not
% affect the memory that is reported.
%
%   memtic
%   n = 125; x = cell(1,n);
%   for i=1:n
%     x{i} = randn(1000,1000); % 8kB per item
%     disp(i);
%   end
%   whos x
%   memtoc
%
% See also TIC, TOC

persistent state

if nargin<1
  action = 'tic';
end

% the memtic/memtoc functions make use of a low-level mex file that interacts directly with the operating system
% do not fail if the mex file does not exist
if isempty(strfind(which('memprofile'), mexext))
  switch action
    case 'tic'
      if nargout
        output = nan;
      end
    case 'toc'
      if nargout
        output = nan;
      end
    otherwise
      error('invalid input argument #1');
  end % switch
  return
end % if mex file does not exist

% get the current usage, this will start memprofile if needed
memstat = memprofile('info');

switch action
  case 'tic'
    if nargin>1
      error('the counter cannot be specified as imput argument');
    end
    
    counter = length(state)+1;
    state(counter).running = true;
    state(counter).time    = memstat(end).time;
    state(counter).mem     = memstat(end).mem;
    
    if nargout
      output = counter;
    end
    
  case 'toc'
    if nargin<2
      % take the latest
      counter = length(state);
    elseif counter<1 || counter>numel(state)
      error('invalid counter');
    end
    
    if counter==0
      memused = 0;
    else
      sel = [memstat.time]>=state(counter).time;
      memused = max([memstat(sel).mem]);       % determine the maximum
      memused = memused - state(counter).mem;  % subtract the initial value
    end
    
    if nargout
      output = memused;
    end
    
  otherwise
    error('invalid input argument #1');
end % switch
