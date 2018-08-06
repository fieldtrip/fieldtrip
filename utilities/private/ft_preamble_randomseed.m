% FT_PREAMBLE_RANDOMSEED is a helper script that records, and optionally sets, the
% random seed/state of the random number generators in RAND/RANDN/RANDI. It calls the
% FieldTrip private RANDOMSEED function, which deals with the MATLAB version
% dependencies.
%
% If cfg.randomseed does not exist, it is set to empty, which indicates default
% behavior.
%
% Default behavior is to only record the current state, not set it (except for MATLAB
% versions prior to 7.3 which set it to a new random state).
%
% Use as
%   ft_preamble randomseed
%   .... regular code goes here ...
%   ft_postamble randomseed
%
% See also FT_POSTAMBLE_RANDOMSEED

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

cfg.randomseed   = ft_getopt(cfg, 'randomseed', []);
ftFuncRandomseed = randomseed(cfg.randomseed);
