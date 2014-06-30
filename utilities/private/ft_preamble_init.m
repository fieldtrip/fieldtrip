% FT_PREAMBLE_INIT is a helper script that display the calling-function's
% help in case the user did not specify any input argument. This can be
% used in all fieldtrip main functions that take at least a cfg input
% argument, and most also take one or multiple data structures.
%
% Use as
%   ft_preamble init

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if nargin==0
  stack = dbstack('-completenames');
  % stack(1) is this script
  % stack(2) is the calling ft_postamble function
  % stack(3) is the main FieldTrip function that we are interested in
  stack = stack(3);
  help(stack.name);
  % throw the error as if it happened in the original function
  msg.message     = 'This function requires one or multiple input arguments, please refer to the documentation above';
  msg.identifier  = '';
  msg.stack       = stack;
  error(msg);
end

% determine whether function execution should be aborted or continued
if isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile) && isfield(cfg, 'outputfilepresent') && ~isempty(cfg.outputfilepresent)
  % check whether the output file already exists
  if ~exist(cfg.outputfile, 'file')
    abort = false;
  else
    % the output file exists, determine how to deal with it
    switch cfg.outputfilepresent
      case 'keep'
        if nargout>0
          % continue executing the parent function
          warning('output file %s is already present, but you also requested an output argument: continuing function execution', cfg.outputfile);
          abort = false;
        else
          % stop executing the parent function
          warning('output file %s is already present: aborting function execution', cfg.outputfile);
          abort = true;
        end
      case 'error'
        error('output file %s is already present', cfg.outputfile);
      case 'overwrite'
        warning('output file %s is already present: it will be overwritten', cfg.outputfile);
        abort = false;
      otherwise
        error('invalid option for cfg.outputfilepresent');
    end % case
  end
else
  abort = false;
end

