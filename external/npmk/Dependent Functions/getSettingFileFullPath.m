function settingFileFullPath = getSettingFileFullPath(functionName)

settingFileFullPath = which(functionName);

if ~isempty(settingFileFullPath)
    settingFilePath = fileparts(settingFileFullPath);
    if ismac
        settingFileFullPath = [settingFilePath '/' functionName '.ini'];
    else
        settingFileFullPath = [settingFilePath '\' functionName '.ini'];
    end
else
    disp('Function not valid.');
    settingFileFullPath = [];
end