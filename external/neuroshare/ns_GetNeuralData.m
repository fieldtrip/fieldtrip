function [ns_RESULT, Data] = ns_GetNeuralData(hFile, EntityID, StartIndex, IndexCount);

%ns_GetNeuralData   Retrieves neural event data by index
%
%   Usage:
%      [ns_RESULT, Data] = 
%               ns_GetNeuralData(hFile, EntityID, StartIndex, IndexCount)
%
%   Description:
%       Returns an array of timestamps for the neural events of the entity
%       specified by EntityID and referenced by the file handle hFile.
%       The index of the first timestamp is StartIndex and the requested
%       number of timestamps is given by IndexCount.  The timestamps are
%       returned in Data.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID	Identification number of the entity in the data file.
%       StartIndex	First index number of the requested Neural Events
%                   timestamp.
%       IndexCount	Number of timestamps to retrieve.
%
%   Return Values:
%       Data	    Array of double precision timestamps.
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

[ns_RESULT, Data] = mexprog(13, hFile, EntityID - 1, StartIndex - 1, IndexCount);
