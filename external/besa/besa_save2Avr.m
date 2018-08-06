function status = besa_save2Avr(file_path, file_name, data_matrix, ...
    time_samples, channel_labels, data_scale_factor, time_scale_factor)
% BESA_SAVE2AVR writes a data matrix into ASCII-vectorized file format 
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
% 
%     [data_matrix]
%         A 2D data matrix to be saved. It should have the size
%         [NumberOfChannels x NumberOfTimeSamples].
% 
%     [time_samples]
%         An array with the time samples corresponding to the data. The
%         size of the array should be [1 x NumberOfTimeSamples].
%
%     [channel_labels]
%         An array containing the channel labels for the current data.    
% 
%     [data_scale_factor]
%         A constant used for converting the data values from Tesla into fT
%         (10e15) or from Volt into uV (10e6). If the data is already in fT
%         or uV the just set it to 1.
% 
%     [time_scale_factor]
%         A constant used for converting the time samples from seconds into
%         milliseconds (1000). If the time samples are already in
%         milliseconds just set the constant to 1.
% 
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 

% Copyright (C) 2013, BESA GmbH
%
% File name: besa_save2Avr.m
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
% Created: 2012-04-17

n_tbins = size(time_samples, 2);

% Check if the matrix is [NumberOfChannels x NumberOfTimeSamples]
% or [NumberOfTimeSamples x NumberOfChannels] instead and if required
% transpose the matrix
tmp_size = size(data_matrix, 1);
% Case [NumberOfTimeSamples x NumberOfChannels]
if(n_tbins == tmp_size)
    
    % Transpose
    data_matrix = data_matrix';
    
end

n_channels = size(data_matrix, 1);

% Scale the time.
start_time = time_scale_factor * time_samples(1);
time_step = time_scale_factor * (time_samples(2) - time_samples(1));

% Scale the data.
data_matrix = data_matrix .* data_scale_factor;

% Change to the directory where the data should be saved.
cd(file_path)

% Open file for writing.
fid = fopen(file_name, 'w');

% Write the first line of the header.
fprintf(fid, ...
    'Npts= %d   TSB= %f DI= %f SB= 1.000 SC= 200.0 Nchan= %d\n', ...
    n_tbins, start_time, time_step, n_channels);

if(iscell(channel_labels))
    
    tmp = channel_labels{1};
    for iCh = 2:n_channels
        
        tmp = [tmp ' ' channel_labels{iCh}];
        
    end
    
    channel_labels = tmp;
    
end
% Write the second line of the header.
fprintf(fid, '%s\n', channel_labels);

% Write the data matrix to the file.
for i=1:n_channels
    
    for j=1:n_tbins

        fprintf(fid, '%f ', data_matrix(i, j));

    end
    
    fprintf(fid, '\n');
    
end

fclose(fid);

status = 1;
