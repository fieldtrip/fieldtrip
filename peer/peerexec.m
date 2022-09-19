function [argout, optout] = peerexec(argin, optin)

% PEEREXEC is the low-level function that executes the job on the
% worker. It also tries to change the path and pwd to those on the
% controller and it catches and deals with any errors in the code that
% is executed.
%
% This function should not be called directly.
%
% See also PEERWORKER, PEERCONTROLLER, PEERCELLFUN

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
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
%
% $Id$
% -----------------------------------------------------------------------

% the fexec function is shared between peer and qsub
[argout, optout] = fexec(argin, optin);

