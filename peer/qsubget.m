function varargout = qsubget(jobid, varargin)

% QSUBGET get the output arguments after the remote job has been executed
% using the Torque or SGE batch queue system.
%
% Use as
%   jobid  = qsubfeval(fname, arg1, arg2, ...)
%   argout = qsubget(jobid, ...)
%
% Optional arguments can be specified in key-value pairs and can include
%   StopOnError    = boolean (default = true)
%   timeout        = number, in seconds (default = 1)
%   sleep          = number, in seconds (default = 0.01)
%   output         = string, 'varargout' or 'cell' (default = 'varargout')
%   diary          = string, can be 'always', 'warning', 'error' (default = 'error')
%
% See also QSUBEVAL, QSUBCELLFUN

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

% instead of reimplementing the whole error handling, just use peerget to do the actual work
varargout = peerget(jobid, varargin);