function varargout = qsubcellfun(fname, varargin)

% QSUBCELLFUN applies a function to each element of a cell-array. The
% function execution is done in parallel using the Torque or SGE batch queue system.
%
% Use as
%   argout = qsubcellfun(fname, x1, x2, ...)
%
% This function has a number of optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   UniformOutput  = boolean (default = false)
%   StopOnError    = boolean (default = true)
%   RetryOnError   = number, number of retries for failed jobs expressed as ratio (default = 0.05)
%   diary          = string, can be 'always', 'never', 'warning', 'error' (default = 'error')
%   timreq         = number, initial estimate for the time required to run a single job (default = 3600)
%   mintimreq      = number, minimum time required to run a single job (default is automatic)
%   memreq         = number, initial estimate for the memory required to run a single job (default = 2*1024^3)
%   minmemreq      = number, minimum memory required to run a single job (default is automatic)
%   order          = string, can be 'random' or 'original' (default = 'random')
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = qsubcellfun(fname, x1, x2);
%
% See also QSUBFEVAL, CELLFUN, FEVAL, DFEVAL, DFEVALASYNC

% -----------------------------------------------------------------------
% Copyright (C) 2011, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

% instead of reimplementing the whole bookkeeping of the jobs, use the peercellfun function to do the actual work
