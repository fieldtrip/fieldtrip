% FT_PREAMBLE_PROVENANCE is a helper script that records the time and memory at the
% start of the function. This is to be used together with FT_POSTAMBLE_PROVENANCE which
% will record and store the time and memory at the end of the function. This is
% stored in the output configuration together with information about the enbvironment,
% such as the name of the user and computer, the matlab and fieldtrip version, etc.
%
% Another aspects of provenance relates to uniquely identifying the input and the
% output data. The code that deals with tracking the information about the input data
% structures is found in ft_preamble_loadvar. The code that deals with tracking the
% information about the output data structures is found in ft_preamble_history.

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

% Record the start time and memory. These are used by ft_postamble_callinfo, which
% stores them in the output cfg.callinfo.  In the mean time, they are stored in the
% function workspace, which is why they should have cryptical names to prevent any
% variable name clashes.

if isfield(cfg, 'trackcallinfo') && ~istrue(cfg.trackcallinfo)
  % do not track the call information
  return
end

ftohDiW7th_FuncTimer = tic();
ftohDiW7th_FuncClock = clock();
ftohDiW7th_FuncMem   = memtic();

