% mff_fileio_read_data() - import MFF data
%
% Usage:
%   >> dat = mff_fileio_read_header(filename);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'header'    - FILEIO structure header
%
% Outputs:
%   dat    - data over the specified range

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function dat = mff_fileio_read_data(filename, varargin)

if nargin < 1
    help mff_fileio_read_data;
    return;
end

header    = ft_getopt(varargin, 'header');

if isempty(header)
    header = mff_fileio_read_header(filename);
end

dat = header.orig.data;

