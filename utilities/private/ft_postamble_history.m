% FT_POSTAMBLE_HISTORY stores the configuration structure that is present in the
% present workspace in the output variable. Furthermore, this function computes the
% MD5 hash of the output data structures for provenance.
%
% Use as
%   ft_postamble history outputvar
%
% See also FT_POSTAMBLE_PROVENANCE

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

% the following section deals with tracking the information about the output data structures
% the corresponding section for the input data structures is in ft_postamble_loadvar

if isfield(cfg, 'trackdatainfo') && istrue(cfg.trackdatainfo)
  % track the information about the output data structures
  if isequal(ft_default.postamble, {'varargout'})
    for i=1:length(varargout)
      % store the hash for each output argument
      cfg.datainfo.output{i} = hashvar(varargout{i});
    end
  else
    for i=1:length(ft_default.postamble)
      cfg.datainfo.output{i} = eval(sprintf('hashvar(%s)', ft_default.postamble{i}));
    end
  end
end

for tmpindx=1:length(ft_default.postamble)
  if isequal(ft_default.postamble, {'varargout'})
    eval(sprintf('try, %s{%d}.cfg = cfg; end', ft_default.postamble{tmpindx}, tmpindx));
  else
    eval(sprintf('try, %s.cfg = cfg; end', ft_default.postamble{tmpindx}));
  end
end
clear tmpindx

% clear warnings from ft_default, so that they don't end up in the next cfg
ft_warning('-clear');
