%% mff_getObject.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%  
%  Returns objects of various types. Some objects are accessed through the
%  'old school' method using javaObject and unmarshall. Others are accessed
%  through the newer method, using factory.openResourceAtURI. 
%%
function theObject = mff_getObject(objType, filename, path)
URI = path;
if objType ~= com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile
    URI = [URI filesep filename];
end
    
oldSchool = true;
switch objType
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info
        objStr = 'com.egi.services.mff.api.Info';
    case com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout
        objStr = 'com.egi.services.mff.api.SensorLayout';
    otherwise
        oldSchool = false;
end
if oldSchool
    theObject = javaObject(objStr, true);
    theObject = theObject.unmarshal(URI, true);
else
    delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
    factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
    resourceVal = objType;
    resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', resourceVal);
%    fprintf('%s %s\n', char(URI), char(resourceType));
    theObject = factory.openResourceAtURI(URI, resourceType);
    if ~isempty(theObject)
        theObject.loadResource();
    end
end

% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info = 7
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack = 3
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories = 9
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal = 2
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs = 4
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN = 8
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout = 10
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Coordinates = 11
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_History = 6
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile = 1
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Photogrammetry = 12
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject = 5
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Any = 0
% com.egi.services.mff.api.MFFResourceType.kMFF_RT_Unknown = -1
