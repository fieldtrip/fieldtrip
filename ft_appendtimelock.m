function [timelock] = ft_appendtimelock(cfg, varargin)

% FT_APPENDTIMELOCK concatenates multiple timelock (ERP/ERF) data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendtimelock(cfg, timelock1, timelock2, ...)
%
% See also FT_TIMELOCKANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble defaults
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% use a helper function to select the consistent parts of the data and to concatenate it
timelock = ft_selectdata(varargin{:}, 'param', {'avg' 'trial' 'cov' 'var' 'dof'});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous varargin
ft_postamble history timelock
ft_postamble savevar timelock
