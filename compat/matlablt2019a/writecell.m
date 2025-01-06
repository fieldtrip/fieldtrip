function writecell(cellarray, filename, varargin)

% This is a very limited overloaded version of the MathWorks function WRITECELL
% which has been introduced in MATLAB 2019a. This version does not offer full
% backwards compatibility, it is only capable of writing text files using
% the (now deprecated) MathWorks function DLMWRITE under the hood.
%
% Use as:
%   writecell(cellarray, filename, varargin)
%
% Please read the code to see which 'varargin' arguments are functional

% Copyright (C) 2022 Jan-Mathijs Schoffelen, Robert Oostenveld
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

sel = find(strcmpi(varargin, 'FileType'));
if ~isempty(sel)
  % remove the new-style argument
  varargin([sel sel+1]) = [];
end

sel = find(strcmpi(varargin, 'WriteMode'));
if ~isempty(sel)
  % remove the new-style argument
  varargin([sel sel+1]) = [];
  % prepend the old-style argument
  varargin = ['-append', varargin];
end

dlmwrite(filename, cellarray, varargin{:})
