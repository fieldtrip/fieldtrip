function status = besa_save2ChanLabels(file_path, file_name, ...
    channel_labels, channel_types)
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
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_save2ChanLabels.m
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

n_channels = n_channelLabels;
if n_channelLabels ~= n_channelTypes
    n_channels = min(n_channelLabels, n_channelTypes);
end;

% Change to the directory where the data should be saved.
cd(file_path)

% Open file for writing.
fid = fopen(file_name, 'w');

% Write the channel labels and the types to the file.
for i=1:n_channels
    switch channel_types{i}
        case 'megmag'
            curChanType = 'MEG';
        case 'eeg'
            curChanType = 'EEG';
        otherwise
            curChanType = 'PGR';
    end

    fprintf(fid, '%s %s\n', curChanType, channel_labels{i});

end

fclose(fid);
status = 1;
