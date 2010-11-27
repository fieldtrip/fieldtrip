function [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID);

%ns_GetAnalogInfo   Retrieves information specific to analog entities
%
%   Usage:
%      [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID)
%
%   Description:
%       Returns information about the Analog Entity associated with
%       EntityID and the file hFile.  The information is stored in 
%       nsAnalogInfo structure.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID		Identification number of the entity in the data
%                       file.
%
%   Return Values:
%       nsAnalogInfo	ns_ANALOGINFO structure to receive the Analog Entity
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

[ns_RESULT, nsAnalogInfo] = mexprog(7, hFile, EntityID - 1);
