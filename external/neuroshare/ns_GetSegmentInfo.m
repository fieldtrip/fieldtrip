function [ns_RESULT, nsSegmentInfo] = ns_GetSegmentInfo(hFile, EntityID);

%ns_GetSegmentInfo   Retrieves information specific to segment entities
%
%   Usage:
%      [ns_RESULT, nsSegmentInfo] = ns_GetSegmentInfo(hFile, EntityID)
%
%   Description:
%       Retrieves information on the Segment Entity, EntityID, in the
%       file referenced by the handle hFile. The information is written to
%       nsSegmentInfo. 
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID		Identification number of the entity in the data
%                       file.
%
%   Return Values:
%       nsSegmentInfo	ns_SEGMENTINFO structure that receives segment
%                       information for the requested Segment Entity.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADFILE	    Invalid file handle passed to 
%                                       function
%                       ns_BADENTITY	Invalid or inappropriate entity 
%                                       identifier specified
%                       ns_FILEERROR	File access or read error
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 6/20/2003

[ns_RESULT, nsSegmentInfo] = mexprog(9, hFile, EntityID - 1);
