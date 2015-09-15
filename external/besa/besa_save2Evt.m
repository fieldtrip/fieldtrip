function status = besa_save2Evt(file_path, file_name, rawdata)

% BESA_SAVE2AEVT writes events into BESA-friendly ASCII file format 
%
% Parameters:
%     [file_path]
%         Full path to the folder where the file should be saved.
% 
%     [file_name]
%         The name of the file where the output should be written.
%
%     [rawdata]
%         General data structure from fieldtrip. Necessary if triggers
%		  with equal spacing should be created. 
% Return:
%     [status] 
%         The status of the writing process. Still unused.
% 
% Author: Robert Spangler
% Created: 2013-10-14

% Change to the directory where the data should be saved.
cd(file_path)

% General information from the rawdata struct
% Number of trials
n_trials = length(rawdata.trial);
% Sampling rate
n_samplingrate = rawdata.fsample;
% Number of samples per trial
n_samplespertrial = size(rawdata.trial{1}, 2);
% Scaling factors to convert time samples from secs to microseconds
time_scale_factor = 1000000;

% Open file for writing.
fid = fopen(file_name, 'w');

% Write first line in file
fprintf(fid, 'Tmu         \tCode\tTriNo\tComnt\n');

% Write all new segment events
for cnt_trial = 1:n_trials
    % starttime of current trial in microseconds
    Tms = (n_samplespertrial*(cnt_trial-1))/n_samplingrate*time_scale_factor;
    % Trigger code
    % 1 = Trigger; 2 = Comment; 3 = Marker; 11 = Pattern1; 12 = Pattern2;
    % 13 = Pattern3; 14 = Pattern4; 15 = Pattern5; 21 = Artifact on;
    % 22 = Artifact off; 31 = Epoch on; 32 = Epoch off; 41 = New segment;
    % 42 = Average segment;   
    Code = '1';
    % Trigger number
    % Trigger: trigger value
    % New segment: date and time of current trial in
    %              YYYY-MM-DDThh-mm-ss.msms [e.g.: 2013-10-14T16:41:00.00]
    TriNo = '1';
    %DateTime = clock;
    %TriNo = sprintf('%i-%i-%iT%02i:%02i:%02.2f',...
    %    DateTime(1), DateTime(2), DateTime(3), ...
    %    DateTime(4), DateTime(5), DateTime(2));
    % comment for current trial
    Comment = sprintf('T#%s', TriNo);
    
    % Write current new segment event to file
    fprintf(fid, '%.0f        	%s\t%s\t%s\n', Tms, Code, TriNo, Comment);    
end

% Close file
fclose(fid);
status = 1;




