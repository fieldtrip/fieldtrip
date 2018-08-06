% FT_POSTAMBLE_PREVIOUS adds the cfg structure from the input data
% structure or structures to the current configuration structure.
%
% Use as
%   ft_postamble previous inputvar
%   ft_postamble previous inputvar1 inputvar2
%   ft_postamble previous varargin

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

% remember the cfg history of the input data structures
cfg.previous = {};
cnt = 0;
for tmpindx=1:length(ft_default.postamble)
  if exist(ft_default.postamble{tmpindx}, 'var')
    tmpvar = eval(ft_default.postamble{tmpindx});
  else
    tmpvar = [];
  end
  if isa(tmpvar, 'struct')
    % the variable is a data structure
    cnt=cnt+1;
    if isfield(tmpvar, 'cfg')
      cfg.previous{cnt} = tmpvar.cfg;
    else
      cfg.previous{cnt} = [];
    end
  elseif isa(tmpvar, 'cell')
    % the variable is a cell-array (i.e. most likely a varargin)
    for cellindx=1:numel(tmpvar)
      cnt=cnt+1;
      if isfield(tmpvar{cellindx}, 'cfg')
        cfg.previous{cnt} = tmpvar{cellindx}.cfg;
      else
        cfg.previous{cnt} = [];
      end
    end
  end
end
clear tmpindx tmpvar cellindx cnt

if length(cfg.previous)==1
  % replace the cell-array by the single struct
  cfg.previous = cfg.previous{1};
end
