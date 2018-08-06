function status = besa_save2Dataset(FullFileName, DataFile, ElpFile, ...
    ControlFile, NumChannels, NumTimeSamples, SamplingInterval, ...
    FirstTimeSample)

% BESA_SAVE2DATASET writes a dataset file for BESA Plot. The file has to 
% have the same name as the DataFile with the  extension 'dataset'. For
% example if the data file is called 'test1.avr' then the dataset file
% should be called 'test1.avr.dataset'.
%
% Parameters:
%     [FullFileName]
%         Full path inclusive file name of the file which is going to be 
%         created.
% 
%     [DataFile]
%         A string containing the name of the data file (e.g. 'test1.avr').
% 
%     [ElpFile]
%         A string containing the name of the file which contains the 
%         channel labels and coordinates (e.g. 'test1.elp').
% 
%     [ControlFile]
%         A string containing the name of the control file for BESA Plot 
%         (e.g. 'test1.bpctrl').   
% 
%     [NumChannels]
%         Number of channels in the data file.
% 
%     [NumTimeSamples]
%         Number of time samples in the data file.
% 
%     [SamplingInterval]
%         Time sampling interval in milliseconds of the data in the data 
%         file.
% 
%     [FirstTimeSample]
%         The time in milliseconds for the first time sample in the data
%         file.
% 
% 
% Return:
%     [status] 
%         The status of the writing process: 1 if the process was 
%         successful and less than 1 if not.
% 

% Copyright (C) 2016, BESA GmbH
%
% File name: besa_save2Dataset.m
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
% Created: 2016-07-06

status = 1;

% Open file for writing.
fid = fopen(FullFileName, 'w');

% MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
% standard output (the screen), and standard error, respectively. When 
% fopen successfully opens a file, it returns a file identifier greater 
% than or equal to 3.
if(fid >= 3)

    fprintf(fid, '[file]\n');
    fprintf(fid, 'dataname=./%s\n', DataFile);
    fprintf(fid, 'sensorname=./%s\n', ElpFile);
    fprintf(fid, 'type=avr\n');
    fprintf(fid, 'conditions=1\n');
    fprintf(fid, '\n');
    fprintf(fid, '[recent]\n');
    fprintf(fid, 'recentControlFileList=./%s\n', ControlFile);
    fprintf(fid, '\n');
    fprintf(fid, '[parameters]\n');
    fprintf(fid, 'blocks=1\n');
    fprintf(fid, 'rows=%i\n', NumChannels);
    fprintf(fid, 'columns=%i\n', NumTimeSamples);
    fprintf(fid, 'channels=%i\n', NumChannels);
    fprintf(fid, 'Samples=%i\n', NumTimeSamples);
    fprintf(fid, 'samplingInterval=%f\n', SamplingInterval);
    fprintf(fid, 'zerotime=%f\n', FirstTimeSample);

    fclose(fid);

else

    status = 0;
    disp('Error! Invalid dataset file identifier.')

end