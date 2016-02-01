%% mff_getObject.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 12/3/2013
%  Copyright 2012, 2013 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%  
%  Returns objects of various types: 
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Any
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Coordinates
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_History
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Photogrammetry
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_PNSSet
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject
%    com.egi.services.mff.api.MFFResourceType.kMFF_RT_Unknown
%%
function theObject = mff_getObject(objType, filename, path)
URI = path;
if objType ~= com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile
    URI = [URI filesep filename];
end
    
delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
resourceVal = objType;
resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', resourceVal);
%    fprintf('%s %s\n', char(URI), char(resourceType));
theObject = factory.openResourceAtURI(URI, resourceType);
if ~isempty(theObject)
    theObject.loadResource();
end
    
