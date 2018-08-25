% FT_POSTAMBLE_SAVEVAR is a helper script that optionally saves the output
% FieldTrip data structures to a *.mat file on disk. This is useful for
% batching and for distributed processing. This makes use of the
% cfg.outputfile variable.
%
% Use as
%   ft_postamble savevar data
%   ft_postamble savevar source mri
%
% See also FT_PREAMBLE_LOADVAR

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
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

global ft_default

% the output data should be saved to a MATLAB file
if isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end

  if isequal(ft_default.postamble, {'varargout'}) && ~iscell(cfg.outputfile)
    % this should be a cell-array, oterwise it cannot be matched with varargout
    cfg.outputfile = {cfg.outputfile};
  end

  if iscell(cfg.outputfile)
    if isequal(ft_default.postamble, {'varargout'})
      % the output is in varargout
      for i=1:length(cfg.outputfile)
        savevar(cfg.outputfile{i}, 'data', varargout{i});
      end % for
      clear i
    else
      % ft_default.postamble is a cell-array containing the variable names
      for i=1:length(cfg.outputfile)
        savevar(cfg.outputfile, ft_default.postamble{i}, eval(ft_default.postamble{i}));
      end % for
      clear i
    end
  else
    % ft_default.postamble{1} contains the name of the only variable
    savevar(cfg.outputfile, ft_default.postamble{1}, eval(ft_default.postamble{1}));
  end
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
  if ~ft_nargout
    % do not return the output variable "ans"
    clear(ft_default.postamble{1});
  end
end

