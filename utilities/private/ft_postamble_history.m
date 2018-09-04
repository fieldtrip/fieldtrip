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

% some fields are for internal use only and should not be stored
cfg = removefields(cfg, ignorefields('history'));

if isequal(iW1aenge_postamble, {'varargout'}) && isequal(iW1aenge_preamble, {'varargin'}) && isfield(cfg, 'previous')
  % distribute the elements of cfg.previous over the output variables
  if iscell(cfg.previous)
    aa5mo0Ke = cfg.previous;
  else
    % this happens when varargin only has a single element
    aa5mo0Ke = {cfg.previous};
  end
  for tmpindx=1:numel(varargout)
    cfg.previous = aa5mo0Ke{tmpindx};
    eval(sprintf('try, varargout{%d}.cfg = cfg; end', tmpindx));
  end
  cfg.previous = aa5mo0Ke;
  clear aa5mo0Ke tmpindx
elseif isequal(iW1aenge_postamble, {'varargout'})
  for tmpindx=1:numel(varargout)
    eval(sprintf('try, varargout{%d}.cfg = cfg; end', tmpindx));
  end
  clear tmpindx
else
  for tmpindx=1:length(iW1aenge_postamble)
    eval(sprintf('try, %s.cfg = cfg; end', iW1aenge_postamble{tmpindx}));
  end
  clear tmpindx
end

