function status = besa_save2PmgPos(file_path, file_name, ...
    channel_coords, channel_orientations, channel_labels)
% besa_save2PmgPos writes the channel coordinates and the corresponding
% labels into the file 'file_name'.
% Brief *.pos and *.pmg file format description:
% - one sensor per line
% - magnetometers: label (optional), six coordinates per line (location,
%   orientation)
% - gradiometers: label (optional), nine coordinates per line (location of
%   primary sensor, location of secondary sensor, orientation).
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
% 
%     [channel_coords]
%         Cell containing the channel coordinates.
% 
%     [channel_orientations]
%         Cell containing the channel coordinates.
% 
%     [channel_labels]
%         Optional parameter. Cell (size: number of channels x 1)
%         containing the channel labels.
%
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_save2PmgPos.m
%
% This file is part of MATLAB2BESA.
%
%    MATLAB2BESA is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    MATLAB2BESA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with MATLAB2BESA. If not, see <http://www.gnu.org/licenses/>.
%
% Author: Robert Spangler
% Created: 2014-01-11

n_channelCoords = size(channel_coords, 1);
n_channelOrientations = size(channel_orientations, 1);
% Check if number of coordinates and orientations match
if n_channelCoords ~= n_channelOrientations
    error('Mismatch in number of sensor coordinates and orientations!');
end;

% Check number of orientations is 3
if size(channel_orientations, 2) ~= 3
    error('Three orientations per sensor required!');
end;

% Check if number of coordinates and labels match
n_channelLabels = 0;
if exist('channel_labels', 'var')
    n_channelLabels = size(channel_labels, 1);
    % Check if number of channel coordinates matches number of channel
    % labels. If they do not match, do not write any labels to file.
    if n_channelCoords ~= n_channelLabels
        n_channelLabels = 0;
    end;
end;

% Change to the directory where the data should be saved.
cd(file_path)

% Open file for writing.
fid = fopen(file_name, 'w');

% Write the channel coordinates and the labels to the file.
for i=1:n_channelCoords
    % Write label
    if n_channelLabels > 0
        fprintf(fid, '%s\t', channel_labels{i});
    end;
    % Write coordinates
    if size(channel_coords{i}, 2) == 3
        % Write magnetometer info: three coordinates
        fprintf(fid, '%2.6f\t%2.6f\t%2.6f\t', channel_coords{i});
    elseif size(channel_coords{i}, 2) == 6
        % Write gradiometer info: six coordinates
        fprintf(fid, '%2.6f\t%2.6f\t%2.6f\t%2.6f\t%2.6f\t%2.6f\t',...
            channel_coords{i});
    else
        % Number of coordinates and orientations does not mach file format
        % description!
        error('Invalid number of coordinates and/or orientations!');
    end;
    
     % Write orientations
     fprintf(fid, '%2.6f\t%2.6f\t%2.6f\n', channel_orientations(i, :));
end

fclose(fid);
status = 1;
