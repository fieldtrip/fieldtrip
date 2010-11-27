function [ns_RESULT, Time] = ns_GetTimeByIndex(hFile, EntityID, Index);

%ns_GetTimeByIndex   Retrieves time range from entity indexes
%
%   Usage:
%      [ns_RESULT, Time] = ns_GetTimeByIndex(hFile, EntityID, Index)
%
%   Description:
%       Retrieves the timestamp for the entity identified by EntityID and
%       numbered Index, from the data file referenced by hFile. The
%       timestamp is returned in Time.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID	Identification number of the entity in the data file.
%       Index		Index of the requested data.
%
%   Return Values:
%       Time	    Variable to receive the timestamp.
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
%   Last modification: 6/20/2003

[ns_RESULT, Time] = mexprog(16, hFile, EntityID - 1, Index - 1);
