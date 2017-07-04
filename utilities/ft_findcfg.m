function [val, status] = ft_findcfg(cfg, var)

% FT_FINDCFG searches for an element in the cfg structure
% or in the nested previous cfgs
%
% Use as
%   val = ft_findcfg(cfg, var)
% where the name of the variable should be specified as string.
%
% e.g.
%   trl   = ft_findcfg(cfg, 'trl')
%   event = ft_findcfg(cfg, 'event')
%
% See also FT_GETOPT, FT_CFG2KEYVAL

% Copyright (C) 2006, Robert Oostenveld
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

% if var(1)~='.'
%   var = ['.' var];
% end
val   = [];
depth = 0;
status = 0;

while ~status
  depth = depth + 1;
  if ~isempty(cfg)
    if issubfield(cfg,  var)
      val = getsubfield(cfg, var);
      status = 1;
    elseif issubfield(cfg, 'previous');
      [val, status] = ft_findcfg(cfg.previous, var);
      if status, break; end
    elseif iscell(cfg)
      for i=1:length(cfg)
        [val, status] = ft_findcfg(cfg{i}, var);
        if status, break; end
      end
    else
      status = -1;
      break
    end
  else
    status = -1;
    break
  end
end
