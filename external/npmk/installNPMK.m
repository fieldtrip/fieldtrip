%
% This function installs Blackrock Microsystems Neural Processing MATLAB
% Kit in MATLAB. Run "help FUNCTIONAME" to get more information about each
% function.
%
% Kian Torab
% kian@blackrockmicro.com
% Blackrock Microsystems
% Version 1.0.0.0
%

function installNPMK

disp('Removing the previous version of Neural Processing MATLAB Kit (NPMK)...');

splitPath = regexp(path, ':', 'split');

for folderIDX = 1:size(splitPath, 2)
    if ~isempty(findstr(splitPath{folderIDX}, 'NPMK'))
        rmpath(splitPath{folderIDX});
    end    
end

disp('Uninstall complete.'); pause(1);

disp('Installing Neural Processing MATLAB Kit (NPMK)...');

currentPath = what;
folderNames = dir;
folderNames(1:2) = [];

addpath(currentPath.path);

for folderIDX = 1:size(folderNames, 1)
    addpath(currentPath.path);
    if folderNames(folderIDX).isdir && folderNames(folderIDX).name(1) ~= '@' && folderNames(folderIDX).name(1) ~= '.'
        addpath(fullfile(currentPath.path, folderNames(folderIDX).name));
    end    
end

savepath;

disp('Install is complete. Run "help FUNCTIONAME" to see how to use each function.');