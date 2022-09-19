function [dataFilename, dataFolder] = getFile(filterSpec, title)

% getFile
%
% Similar to uigetfile, it opens a dialogue to get a file, but it also
% remembers the last folder location, so the user will not need to navigate
% to the same folder multiple times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use [dataFilename dataFolder] = getFile(filterSpec)
%
%   filterSpec:   The type of file to show. For example, if set to '*.txt'
%   (OPTIONAL)    only text files will be visible for choosing.
%                 DEFAULT: *.*, so all files are shown.
%
%   title:        This parameter will be the title of the dialogue box.
%   (OPTIONAL)    DEFAULT: No title will be shown.
%
%   dataFilename: The name of the chosen file.
%
%   dataFolder:   The path of the chosen file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USAGE EXAMPLE: 
%   
%   [fileName pathName] = getFile('*.jpg');
%
%   In this example a dialogue will open and only show all jpg files. The
%   name of the chosen file and the path to it are stores in fileName and
%   pathName respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Salt Lake City, UT
%   Contributors: 
%   
%   Version 1.1.0.0
%   Last edit by: Kian Torab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.1.0:
%   - Initial release.
%
% 1.0.1.0:
%   - Minor update to the help and some formatting. 
%
% 1.0.2.0: July 4, 2014
%   - Allows for selection of NSx files only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('filterSpec', 'var')
    filterSpec = '*.*';
end

if ~exist('title', 'var')
    title = '';
end

settingFileFullPath = getSettingFileFullPath('getFile');

%% Opens the getFolder.ini file to see what the last accessed folder was
if exist(settingFileFullPath, 'file') == 2
    settingsFID = fopen(settingFileFullPath, 'r');
    defaultOpenLocation = fscanf(settingsFID, '%200c');
    fclose(settingsFID);
    backDir = cd;
    if exist(defaultOpenLocation, 'dir') == 7
        cd(defaultOpenLocation);
    end
    [dataFilename, dataFolder] = uigetfile(filterSpec, title);
    cd(backDir);
else
    [dataFilename, dataFolder] = uigetfile(filterSpec, title);
end

%% Saves the last opened folder in the getFolder.ini file for later use
if ischar(dataFolder)
    settingsFID = fopen(settingFileFullPath, 'w');
    fprintf(settingsFID, '%s', dataFolder);
    fclose(settingsFID);
end