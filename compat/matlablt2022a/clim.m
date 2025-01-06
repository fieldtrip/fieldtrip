function clim(varargin)

% This is an drop-in replacement version of the MathWorks function CLIM
% which has been introduced in MATLAB 2022a. This replacement function
% allows CLIM also to be called on older MATLAB versions.
%
% The directory containing this function should only be added to the path
% of MATLAB versions prior to 2022a.

% Copyright (C) 2024, Robert Oostenveld
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

caxis(varargin{:});
