% FT_POSTAMBLE_HASTOOLBOX is executed at the end of each FieldTrip
% function to remove other toolboxes that have been added automatically
% by FT_HASTOOLBOX during execution of the specific function.
%
% Use as
%   ft_postamble hastoolbox
%
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_HASTOOLBOX

% Copyright (C) 2018, Robert Oostenveld, DCCN
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

if ~isempty(ft_default) && isfield(ft_default, 'toolbox') && isfield(ft_default.toolbox, 'cleanup')
  while ~isempty(ft_default.toolbox.cleanup)
    toolbox = ft_default.toolbox.cleanup{end};
    ft_warning('off','backtrace');
    ft_warning('removing %s toolbox from your MATLAB path', toolbox);
    ft_warning('on','backtrace');
    rmpath(genpath(toolbox));
    ft_default.toolbox.cleanup = ft_default.toolbox.cleanup(1:end-1);
  end
end
