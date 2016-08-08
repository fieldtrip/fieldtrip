function status = besa_save2Elp(file_path, file_name, ...
    spherical_coordinates, channel_labels, channel_type)
% BESA_SAVE2ELP writes a matrix of spherical coordinates into ASCII file 
% in the ELP-format. 
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
% 
%     [spherical_coordinates]
%         A 2D matrix containing the coordinates. It should have the size
%         [NumberOfChannels x 2(or 3)]. For each row the structure should
%         be [theta, phi, (r)], the radius r is optional and is not needed
%         for the export, since it is assumed that the coordinates are on
%         unit sphere and r is always 1.
% 
%     [channel_labels]
%         A cell-array of strings containing the channel labels for the 
%         current data. It is also accepted if this a regular array 
%         containing the channel labels separated with whitespaces.   
% 
%     [channel_type]
%         A single string containing the type of the channels. Possible
%         values are: MEG, EEG, POL, ICR.
% 
% 
% Return:
%     [status] 
%         The status of the writing process: 1 if the process was 
%         successful and less than 1 if not.
% 

% Copyright (C) 2015, BESA GmbH
%
% File name: besa_save2Elp.m
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
% Author: Todor Jordanov
% Created: 2015-07-29

status = 1;

NumChannels = size(spherical_coordinates, 1);

FullPath = fullfile(file_path, file_name);

% If the labels are not stored in a cell array then create a cell array of
% them.
if(~iscell(channel_labels))
    
    tmp_cell = strsplit(channel_labels);
    channel_labels = tmp_cell;
    
end

% Open file for writing.
fid = fopen(FullPath, 'w');

% MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
% standard output (the screen), and standard error, respectively. When 
% fopen successfully opens a file, it returns a file identifier greater 
% than or equal to 3.
if(fid >= 3)

    % Write the coordinates to the file.
    for i=1:NumChannels

        fprintf(fid, '%s\t%s\t%f\t%f \n', channel_type, ...
            channel_labels{i}, spherical_coordinates(i, 1), ...
            spherical_coordinates(i, 2));

    end

    fclose(fid);

else
    
    status = -1;
    disp('Error! Invalid file identifier.')
    
end