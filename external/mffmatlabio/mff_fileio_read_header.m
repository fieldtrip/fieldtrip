% mff_fileio_read_header() - import MFF header
%
% Usage:
%   >> header = read_eeglabheader(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure

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

function header = mff_fileio_read_header(mffFile)

if nargin < 1
    help mff_fileio_read_header;
    return;
end

% import basic information
tmp = mff_import(mffFile);
curPath = pwd;
p = fileparts(which('ft_read_header'));
cd(fullfile(p, 'private'));
header = read_eeglabheader( tmp ); 
cd(curPath);
header.orig = tmp;