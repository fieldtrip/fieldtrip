% mff_fileio_read_event() - import MFF events
%
% Usage:
%     >> event = mff_fileio_read_event(filename, ...);
%
% Inputs:
%     filename - [string] file name
%
% Optional inputs:
%     'header' - FILEIO structure header
%
% Outputs:
%     event - FILEIO toolbox event structure

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

function event = mff_fileio_read_event(filename, varargin)

if nargin < 1
    help mff_fileio_read_event;
    return;
end

hdr = ft_getopt(varargin, 'header');

if isempty(hdr)
    hdr = mff_fileio_read_header(filename);
end

curPath = pwd;
p = fileparts(which('ft_read_header'));
cd(fullfile(p, 'private'));

% bypass the last part where time-locking events are added
hdr.orig.trials = 1;

event = read_eeglabevent( [], 'header', hdr, varargin{:} ); 

cd(curPath);
