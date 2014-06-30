%% write_mff_history.m
%  Matlab File
%  author Colin Davey
%  date 4/15/2014
%  Copyright 2014 EGI. All rights reserved.
%
%  This function adds an item to the list of history items. It is intended
%  to be run after running code that either a) creates a new mff file based
%  on an existing one, or b) modifies an existing one. The history item can
%  be seen in Net Station through the File Info->History feature. 
%
%  filePath ? The path to the .mff file. 
%
%  entryStruct has the following elements: 
% 
%  name - string that corresponds to the name of the tool specification in
%  NetStation. For example, if you create a segmentation specification for
%  a VTD experiment, and call it "VTD Seg", then the string "VTD Seg" would
%  go here.
%
%  method - The name of the tool (not the tool specification.) In the above
%  example, the string "Segmentation" would go here.
%
%  version - The version of the tool, eg "1.0".
%
%  beginTime - The time the tool started running. The following matlab
%  command creates a time in the proper format:
%    sprintf('%d-%02d-%02dT%02d:%02d:%02.5f',clock);
%
%  endTime - The time the tool finished running. See beginTime, above, for
%  the matlab command that creates a time in the proper format.
% 
% sourceFileList - The list of files that were processed by the tool to
% generate the resulting file. For example, a tool that filters a single
% file to create a filtered version of the data would have just one item in
% this list. A tool that creates a grand average by averaging multiple
% single-subject files would have multiple items in this list. 
% 
% settingList - A cell array of zero or more strings that express the
% settings/parameters for the tool. These are optional, and are intended to
% be read by humans, not machines.
% 
% resultList - A cell array of zero or more strings that express the
% results of the tool. These are optional, and are intended to be read by
% humans, not machines.
%%
function write_mff_history(filePath, entryStruct)
try
    mff_valid(filePath);
catch theException
    throw(theException);
end

%% 1) Read in history resource.
histObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_History, 'history.xml', filePath);

%% 2) Get the list of current entires (it might be null, in which case
% create a new one and set it).
if isempty(histObj)
    delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
    factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
    resourceVal = com.egi.services.mff.api.MFFResourceType.kMFF_RT_History;
    resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', resourceVal);
    factory.createResourceAtURI([filePath filesep 'history.xml'], resourceType);
    histObj = factory.openResourceAtURI([filePath filesep 'history.xml'], resourceType);
    
    entryList = javaObject('java.util.ArrayList');
    histObj.setEntries(entryList);
else
    entryList = histObj.getEntries;
    if isempty(entryList) % Is this feasible? ie history exists, but entryList is empty? 
        entryList = javaObject('java.util.ArrayList');
        histObj.setEntries(entryList);
    end
end

%% 3) Add an Entry object to the list. 
entry = javaObject('com.egi.services.mff.api.Entry');
tool = javaObject('com.egi.services.mff.api.Tool');
tool.setName(entryStruct.name);
tool.setKind('Transformation');
tool.setMethod(entryStruct.method);
tool.setVersion(entryStruct.version);
tool.setBeginTime(entryStruct.beginTime);
tool.setEndTime(entryStruct.endTime);

sourceFileList = javaObject('java.util.ArrayList');
for p=1:size(entryStruct.sourceFileList,2)
    filePathObj = javaObject('com.egi.services.mff.api.FilePath');
    filePathObj.setFilePath(entryStruct.sourceFileList{p});
    sourceFileList.add(filePathObj);
    filePathObj = [];
end
java.lang.Runtime.getRuntime.freeMemory;
tool.setSourceFiles(sourceFileList);

settingList = javaObject('java.util.ArrayList');
for p=1:size(entryStruct.settingList,2)
    settingList.add(entryStruct.settingList{p});
end
tool.setSettings(settingList);

resultList = javaObject('java.util.ArrayList');
for p=1:size(entryStruct.resultList,2)
    resultList.add(entryStruct.resultList{p});
end
tool.setResults(resultList);

entry.setEntry(tool);
entry.setType('tool');

entryList.add(entry);

%% 4) Save the history resource in the standard way all resources are saved.
histObj.setEntries(entryList);
histObj.saveResource;
