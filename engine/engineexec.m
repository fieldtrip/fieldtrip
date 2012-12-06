function [argout, optout] = engineexec(argin, optin)

% ENGINEEXEC is a helper function to execute a command in another engine.
% Normally you should not start this function yourself, but rather
% use ENGINECELLFUN or ENGINEFEVAL.
%
% See also ENGINEPOOL, ENGINECELLFUN, ENGINEFEVAL, ENGINEGET

% -----------------------------------------------------------------------
% Copyright (C) 2012, Robert Oostenveld
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

try
  % the actual work is being done by fexec, which is shared among
  % the qsub, peer and engine toolboxes
  [argout, optout] = fexec(argin, optin);
  
catch err
  % this is to avoid MATLAB from hanging in case fexec fails, since
  % after the job execution we want MATLAB to exit
  disp(err);
  warning('an error was caught');
  
end % try-catch
