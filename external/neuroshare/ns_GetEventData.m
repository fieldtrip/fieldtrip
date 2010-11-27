function [ns_RESULT, TimeStamp, Data, DataSize] = ns_GetEventData(hFile, EntityID, Index);

%ns_GetEventData   Retrieves event data by index
%
%   Usage:
%      [ns_RESULT, TimeStamp, Data, DataSize] = 
%                                   ns_GetEventData(hFile, EntityID, Index)
%
%   Description:
%       Returns the data values from the file referenced by hFile and the
%       Event Entity EntityID.  The Event data entry specified by Index
%       is written to Data and the timestamp of the entry is returned to
%       TimeStamp.  Upon return of the function, the value at DataSize
%       contains the number of bytes actually written to Data.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID	    Identification number of the entity in the data
%                       file.
%       Index	        The index number of the requested Event data item.
%
%   Return Values:
%       TimeStamp	    Variable that receives the timestamp of the Event
%                       data item.
%       Data	        Variable that receives the data for the Event entry.
%                       The format of Event data is specified by the member
%                       EventType in ns_EVENTINFO.
%       DataSize	    Variable that receives the actual number of bytes
%                       of data retrieved in the data buffer.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADFILE	    Invalid file handle passed to 
%                                       function
%                       ns_BADENTITY	Invalid or inappropriate entity 
%                                       identifier specified
%                       ns_BADINDEX	    Invalid entity index specified
%                       ns_FILEERROR	File access or read error
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 8/11/2003

[ns_RESULT, TimeStamp, Data, DataSize] = mexprog(6, hFile, EntityID - 1, Index - 1);
