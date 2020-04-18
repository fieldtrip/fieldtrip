function settingFileFullPath = getSettingFileFullPath(functionName)

settingFileFullPath = which(functionName);

if ~isempty(settingFileFullPath)
    settingFilePath = fileparts(settingFileFullPath);
    settingFileFullPath = [settingFilePath '\' functionName '.ini'];
else
    disp('Function not valid.');
    settingFileFullPath = [];
end