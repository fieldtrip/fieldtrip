function [ns_RESULT, nsEventInfo] = ns_GetEventInfo(hFile, EntityID);

%ns_GetEventInfo   Retrieves information specific to event entities
%
%   Usage:
%      [ns_RESULT, nsEventInfo] = ns_GetEventInfo(hFile, EntityID)
%
%   Description:
%       Retrieves information from the file referenced by hFile about the
%       Event Entity, EntityID, in the structure nsEventInfo.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID	    Identification number of the entity in the data
%                       file.
%
%   Return Values:
%       nsEventInfo	    ns_EVENTINFO structure to receive the Event Entity
%                       information.
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

[ns_RESULT, nsEventInfo] = mexprog(5, hFile, EntityID - 1);
