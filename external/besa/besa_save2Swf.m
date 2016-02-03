function status = besa_save2Swf(file_path, file_name, data_matrix, ...
    time_samples, source_labels, data_scale_factor, time_scale_factor)
	
% BESA_SAVE2SWF writes a data matrix into ASCII-source waveform file format 
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
%     [source_labels]
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
% File name: besa_save2Swf.m
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
% Created: 2012-08-27

n_sources = size(source_labels, 2);%size(data_matrix, 1);

% Check if data is vectorized or multiplex.
isMul = 1;  % boolean variable containing the info if the file is 
            %vectorized or multiplexed.
tmp1 = size(data_matrix);
if(tmp1(1) == n_sources)
    
    n_tbins = tmp1(2);
    isMul = 0;

else
    
    n_tbins = tmp1(1);
    
end

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
    n_tbins, start_time, time_step, n_sources);

% Write the data matrix to the file.
for i=1:n_sources
    
    fprintf(fid, '%s: ', source_labels{i});
    
    for j=1:n_tbins

        if(isMul == 0)
            
            fprintf(fid, '%f ', data_matrix(i, j));
        
        else
            
            fprintf(fid, '%f ', data_matrix(j, i));
            
        end

    end
    
    fprintf(fid, '\n');
    
end

fclose(fid);

status = 1;
