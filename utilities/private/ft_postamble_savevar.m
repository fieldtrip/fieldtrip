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

global ft_default

% the output data should be saved to a MATLAB file
if isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end
  
  % ft_default.postamble{1} contains the name of the variable
  savevar(cfg.outputfile, ft_default.postamble{1}, eval(ft_default.postamble{1}));
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
  if ~nargout
    % do not return the output variable "ans"
    clear(ft_default.postamble{1});
  end
end
