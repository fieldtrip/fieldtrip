% FT_POSTAMBLE_HISTORY is a helper script that stores the configuration structure that
% is present in the present workspace (i.e. the workspace of the calling function) in
% the output variable.
%
% Use as
%   ft_postamble history outputvar
%
% See also FT_POSTAMBLE_PROVENANCE

% Copyright (C) 2011-2016, Robert Oostenveld, DCCN
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

% some fields are for internal use only and should not be stored
cfg = removefields(cfg, ignorefields('history'));

for tmpindx=1:length(ft_default.postamble)
  if isequal(ft_default.postamble, {'varargout'})
    eval(sprintf('try, %s{%d}.cfg = cfg; end', ft_default.postamble{tmpindx}, tmpindx));
  else
    eval(sprintf('try, %s.cfg = cfg; end', ft_default.postamble{tmpindx}));
  end
end
clear tmpindx
