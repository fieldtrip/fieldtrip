function inspect_pull1419


% The code requires the MatNWB toolbox in a version => 
% NeurodataWithoutBorders/matnwb@a1b1266 (19.05.2020).
% 
% It is assumed that the user is familiar with the basics of MatNWB and 
% the user's system is in a state that would also allow loading the data using MatNWB 
% (MatNWB must be in the path, generateCore must have been run, the correct schema must be 
% installed, see NWB schemas).
% 
% The code was tested with the following example datasets that are freely available online:
%
% From https://gui.dandiarchive.org/#/file-browser/folder/5e6eb2b776569eb93f451f8d:
% sub-YutaMouse20_ses-YutaMouse20-140324_behavior+ecephys.nwb
% 
% From https://osf.io/hv7ja/ 
% P9HMH_NOID5.nwb
%
% At its current stage, the code is not feature-complete. For example:
% 
% we have not implemented ft_read_event. NWB:N is a pretty generic dataformat and can 
% contain very diverse types of data. So we were not sure how to programmatically and 
% reliably create an event output that could be used in a trial function. Waveforms are 
% not read in by ft_read_spike
%
% This script illustrates some simple tests we did on the data using the
% example datasets. 
%
% Folder structure
% ----------------
% |CD
%    |matnwb
%       |nwb-schema
%          |core
%             |<Files in core might need to be substituted with files from the
%               schema-version of the NWB file of interest>
%
% calling generateCore() from matnwb will create the following folder in the CD based
% on the nwb-schema in matnwb. This needs to be done only once, but with
% the right schema version for the data. 
%
% |+types
%
% 
%
% Useful links:
%
% For examples of spike data structures in FieldTrip see
% http://www.fieldtriptoolbox.org/tutorial/spike/


%% Conversion Example 1. NWB data, Version 2.1.0
% from https://osf.io/hv7ja/
% About the data:
% https://github.com/rutishauserlab/recogmem-release-NWB

% Make sure no duplicate versions of +types are in the search path
restoredefaultpath()

% Navigate to working directory (+types lives here, too, after generateCore())
cd('D:/NWB/hackathon_may_2020') % change as needed

% Add fieldtrip to search path and initialize
addpath('../fieldtrip') % change as needed
ft_defaults

% Add matnwb to search path
% matnwb release 2.2.4, with schema 2.1.0 (manually replaced core in nwb-schema)
addpath(genpath('../matnwb-searchFor_schema-2.1.0')) 


% (re-)generate core classes for matNWB from it's schema (lands in +types)
generateCore()

% Filename and path
nwbFile = 'D:\NWB\data\NWBData\NWBData\P9HMH_NOID5.nwb'; 
disp(util.getSchemaVersion(nwbFile))

% Load data in nwb format
nwb = nwbRead(nwbFile);
disp(nwb)

% Try to obtain hdr, lfp data and spike data in FieldTrip format
try
    hdr = ft_read_header(nwbFile); % contains no lfp data: throws error
catch ME
    disp('Could not load in hdr information')
    rethrow(ME)
end
try
    dat = ft_read_data(nwbFile); % contains no lfp data: throws error
catch ME
    disp('Could not load in NWB data')
    rethrow(ME)
end
try
    spike = ft_read_spike(nwbFile); % contains spike data: Converts
catch ME
    disp('Could not read in spike data from NWB file.')
    rethrow(ME)
end



%% Conversion Example 2. NWB data, Version 2.0.1
% 2.0.1
% from https://gui.dandiarchive.org/#/file-browser/folder/5e6eb2b776569eb93f451f8d
% About the data:
% http://www.buzsakilab.com/content/PDFs/Senzai2017Neuron.pdf

% Make sure no duplicate versions of +types are in the search path
restoredefaultpath()

% Navigate to working directory (+types lives here, too, after generateCore())
cd('D:/NWB/hackathon_may_2020') % change as needed

% Add fieldtrip to search path and initialize
addpath('../fieldtrip') % change as needed
ft_defaults

% Add matnwb to search path
% matnwb release 2.2.4, with schema 2.0.1 (manually replaced core in nwb-schema)
% addpath(genpath('../matnwb-searchFor_schema-2.0.1')) 
addpath(genpath('../matnwb-2.2.4_schema-2.0.1'))
% (re-)generate core classes for matNWB from it's schema (lands in +types)
generateCore()

% Filename and path
nwbFile = 'D:\NWB\data\sub-YutaMouse20_ses-YutaMouse20-140324_behavior+ecephys.nwb'; 
disp(util.getSchemaVersion(nwbFile))

% Load data in nwb format
nwb = nwbRead(nwbFile);
disp(nwb)

% Try to obtain hdr, lfp data and spike data in FieldTrip format
try
    hdr = ft_read_header(nwbFile);
catch ME
    disp('Could not load in hdr information')
    rethrow(ME)
end
try
    dat = ft_read_data(nwbFile);
catch ME
    disp('Could not load in NWB data')
    rethrow(ME)
end
try
    spike = ft_read_spike(nwbFile);
catch ME
    disp('Could not read in spike data from NWB file.')
    rethrow(ME)
end


% eof

