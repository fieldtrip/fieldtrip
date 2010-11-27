function [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID);

%ns_GetEntityInfo   Retrieves general entity information and type
%
%   Usage:
%      [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID) 
%
%   Description:
%       Retrieves general information about the entity, EntityID, from 
%       the file referenced by the file handle hFile.  The information is
%       passed in the structure nsEntityInfo.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID	Identification number of the entity in the data file.
%                   The total number of entities in the data file is
%                   provided by the member EntityCount in the ns_FILEINFO
%                   structure.
%
%   Return Values:
%       nsEntityInfo	ns_ENTITYINFO structure to receive entity
%                       information
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

[ns_RESULT, nsEntityInfo] = mexprog(4, hFile, EntityID -1);
