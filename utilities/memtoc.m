function memused = memtoc(counter)

% MEMTOC return the memory that was used
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
% Example: measure the memory increase due to allocating a lot of memory
%   memtic
%   n = 125; x = cell(1,n);
%   for i=1:n
%     x{i} = zeros(1000,1000); % 8kB per item
%     disp(i);
%   end
%   whos x
%   memtoc
%
% See also TIC, TOC

% this function itself does not do anything, because it depends 
% on the persistent variable inside memtic
if nargin==0
  memused = memtic('toc');
else
  memused = memtic('toc', counter);
end

if ~nargout
  fprintf('Estimated memory use is %.03f MB\n', memused/(1024^2));
  clear memused
end
