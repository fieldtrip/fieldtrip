function status = besa_save2ChanPositions(file_path, file_name, ...
    channel_labels, channel_types, channel_positions, coil_orientations)
% besa_save2ChanLabels writes the channel labels and the corresponding type
% into the file 'file_name'.
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
%         A cell (size: number of channels x 1) containing the channel
%         types for the current data.
%
%     [channel_positions]
%         A cell (size: number of channels x 1) containing the channel
%         positions for the current data.
%
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_save2ChanPositions.m
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
% Created: 2013-06-25

n_channelLabels = size(channel_labels, 1);
n_channelTypes = size(channel_types, 1);
n_channelPositions = size(channel_positions, 1);
n_coilOrientations = size(coil_orientations, 1);

n_channels = n_channelLabels;
if (n_channelLabels ~= n_channelTypes ||...
    n_channelLabels ~= n_channelPositions ||...
    n_channelTypes ~= n_channelPositions)    
        n_channels = min(n_channelLabels, n_channelTypes,...
            n_channelPositions);
end;

% Change to the directory where the data should be saved.
cd(file_path)

% Open file for writing.
fid = fopen(file_name, 'w');

% Write the channel labels and the positions to the file.
for i=1:n_channels
    switch channel_types{i}
        case 'megmag'
            xPos = channel_positions(i,2);
            yPos = channel_positions(i,1);
            zPos = channel_positions(i,3);
            xOri = coil_orientations(i,2);
            yOri = coil_orientations(i,1);
            zOri = coil_orientations(i,3);
            fprintf(fid, '%s %f %f %f %f %f %f\n', channel_labels{i},...
                xPos, yPos, zPos, xOri, yOri, zOri);
        case 'eeg'
            %x_pos = channel_positions(i, 1);
            %y_pos = channel_positions(i, 2); 
            %z_pos = channel_positions(i, 3);
            %radius = sqrt(x_pos^2 + y_pos^2 + z_pos^2);
            %inclination = acos(z_pos/radius); 
            %azimuth = atan2(y_pos/x_pos);       
        otherwise
    end   
end

fclose(fid);
status = 1;
