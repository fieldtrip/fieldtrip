function [ns_RESULT, nsSegmentSourceInfo] = ns_GetSegmentSourceInfo(hFile, EntityID, SourceID);

%ns_GetSegmentSourceInfo   Retrieves information about the sources that
%   generated the segment data
%
%   Usage:
%       [ns_RESULT, nsSegmentSourceInfo] = 
%               ns_GetSegmentSourceInfo(hFile, EntityID, SourceID)
%
%   Description:
%       Retrieves information about the source entity, SourceID, for the
%       Segment Entity identified by EntityID, from the file referenced
%       by the handle hFile. The information is written to the
%       nsSegmentSourceInfo.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID	    Identification number of the Segment Entity.
%       SourceID	    Identification number of the Segment Entity source.
%
%   Return Values:
%       nsSegmentSourceInfo     ns_SEGSOURCEINFO structure that receives
%                               information about the source.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADFILE	    Invalid file handle passed to 
%                                       function
%                       ns_BADENTITY	Invalid or inappropriate entity 
%                                       identifier specified
%                       ns_BADSOURCE	Invalid source identifier specified
%                       ns_FILEERROR	File access or read error
%
%   Remarks:
%       The value of SourceID is an integer index value ranging from
%       0 to SourceCount-1 (which is returned by the function
%       ns_GetSegmentInfo).
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 6/20/2003

[ns_RESULT, nsSegmentSourceInfo] = mexprog(10, hFile, EntityID - 1, SourceID - 1);
