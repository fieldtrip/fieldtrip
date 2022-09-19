% getFolder
%
% Similar to uigetdir except that it will remember the location of the
% last opened folder.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use plotAverageWaveforms
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Salt Lake City, UT
%   
%   Version 1.0.0.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataFolder = getFolder

settingFileFullPath = getSettingFileFullPath('getFolder');

%% Opens the getFolder.ini file to see what the last accessed folder was
if exist(settingFileFullPath, 'file') == 2
    settingsFID = fopen(settingFileFullPath, 'r');
    defaultOpenLocation = fscanf(settingsFID, '%200c');
    fclose(settingsFID);
else
    defaultOpenLocation = [];
end

%% Gets a folder by opening the last accessed folder as the last location
dataFolder = uigetdir(defaultOpenLocation);

%% Saves the last opened folder in the getFolder.ini file for later use
if ischar(dataFolder)
    settingsFID = fopen(settingFileFullPath, 'w');
    fprintf(settingsFID, '%s', dataFolder);
    fclose(settingsFID);
end