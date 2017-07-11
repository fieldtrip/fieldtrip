function output = memtic(action, counter)

% MEMTIC start a MATLAB memory recorder
%
%   MEMTIC and MEMTOC functions work together to measure memory usage.
%   MEMTIC, by itself, saves the current memory footprint that MEMTOC
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

% Copyright (C) 2011, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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
      ft_error('invalid input argument #1');
  end % switch
  return
end % if mex file does not exist

% get the current usage, this will start memprofile if needed
memstat = memprofile('info');

switch action
  case 'tic'
    if nargin>1
      ft_error('the counter cannot be specified as imput argument');
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
      ft_error('invalid counter');
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
    ft_error('invalid input argument #1');
end % switch

