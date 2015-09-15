function status = besa_save2SurfacePoint(file_path, file_name, ...
    channel_labels, channel_coords, fiducial_labels, fiducial_coords, ...
    headsurfacepoint_coords)
% besa_save2SurfacePoint writes the channel and fiducial coordinates into
% the file 'file_name'.
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
% 
%     [channel_labels]
%         A cell (size: number of channels x 1) containing the channel
%         labels for the current data.
% 
%     [channel_types]
%         An array (size: number of channels x 3) containing the channel
%         coordinates for the current data.
%
%     [fiducial_labels]
%         A cell (size: number of fiducials x 1) containing the fiducial
%         labels for the current data.
% 
%     [fiducial_coords]
%         An array (size: number of fiducials x 3) containing the fiducials
%         coordinates for the current data.
%
%     [headsurfacepoint_coords]
%         Optional. An array (size: number of fiducials x 3) containing
%         additional head surface point coordinates for the current data.
%         Labels for head surface points will be created automatically,
%         starting with Sfh1.
%
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_save2SurfacePoint.m
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
% Created: 2013-12-22

n_FiducialLabels = size(fiducial_labels, 1);
n_FiducialCoords = size(fiducial_coords, 1);
n_ChannelLabels = size(channel_labels, 1);
n_ChannelCoords = size(channel_coords, 1);

% Get number of channels and fiducials
n_NumFiducials = n_FiducialLabels;
if n_FiducialLabels ~= n_FiducialCoords
    n_NumFiducials = min(n_FiducialLabels, n_FiducialCoords);
end;
n_NumChannels = n_ChannelLabels;
if n_ChannelLabels ~= n_ChannelCoords
    n_NumChannels = min(n_ChannelLabels, n_ChannelCoords);
end;

if exist('headsurfacepoint_coords', 'var');
    n_HSPCoords = size(headsurfacepoint_coords, 1);
else
    n_HSPCoords = 0;
end;
    
% Change to the directory where the data should be saved.
cd(file_path)

% Open file for writing.
fid = fopen(file_name, 'w');

% Write the fiducial labels and the coordinates to the file.
for i=1:n_NumFiducials
    fprintf(fid, '%s\t\t\t%2.6f\t%2.6f\t%2.6f\n', fiducial_labels{i}, ...
        fiducial_coords(i,1), fiducial_coords(i,2), fiducial_coords(i,3));
end;

% Write the channel labels and the coordinates to the file.
for i=1:n_NumChannels
    curChannelLabel = strrep(channel_labels{i}, 'Ele_', '');
    fprintf(fid, '%s\t\t\t%2.6f\t%2.6f\t%2.6f\n', curChannelLabel, ...
        channel_coords(i,1), channel_coords(i,2), channel_coords(i,3));
end

% Write additional head surface point labels and coordinates to the file.
for i=1:n_HSPCoords
    curHSPLabel = strcat('Sfh', num2str(i));
    fprintf(fid, '%s\t\t\t%2.6f\t%2.6f\t%2.6f\n', curHSPLabel, ...
        headsurfacepoint_coords(i,1), headsurfacepoint_coords(i,2), ...
        headsurfacepoint_coords(i,3));
end;

fclose(fid);
status = 1;



