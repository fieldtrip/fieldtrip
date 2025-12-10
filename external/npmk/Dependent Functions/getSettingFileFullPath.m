function settingFileFullPath = getSettingFileFullPath(functionName)
% settingFileFullPath
% 
% Returns the full path of the settings file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0:
%   - Updated to support unix file system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settingFileFullPath = which(functionName);

if ~isempty(settingFileFullPath)
    settingFilePath = fileparts(settingFileFullPath);
    if ismac || isunix
        settingFileFullPath = [settingFilePath '/' functionName '.ini'];
    elseif isunix
        settingFileFullPath = [settingFilePath '/' functionName '.ini'];
    else
        settingFileFullPath = [settingFilePath '\' functionName '.ini'];
    end
else
    disp('Function not valid.');
    settingFileFullPath = [];
end