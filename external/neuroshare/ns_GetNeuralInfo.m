function [ns_RESULT, nsNeuralInfo] = ns_GetNeuralInfo(hFile, EntityID);

%ns_GetNeuralInfo   Retrieves information for neural event entities
%
%   Usage:
%      [ns_RESULT, nsNeuralInfo] = ns_GetNeuralInfo(hFile, EntityID)
%
%   Description:
%       Retrieves information on Neural Event entity EntityID from the
%       file referenced by hFile.  The information is returned in the
%       nsNeuralInfo.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID	Identification number of the entity in the data file.
%
%   Return Values:
%       nsNeuralInfo    ns_NEURALINFO structure to receive the Neural
%                       Event information.
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
%   Last modification: 10/24/2003

[ns_RESULT, nsNeuralInfo] = mexprog(12, hFile, EntityID - 1);

SourceEntityID = [nsNeuralInfo.SourceEntityID] + 1;
SourceUnitID = [nsNeuralInfo.SourceUnitID];

ind = find(SourceUnitID == 1);
SourceUnitID(ind) = 255;

ind = find((SourceUnitID > 1) & (SourceUnitID < 255));
SourceUnitID(ind) = log2(SourceUnitID(ind));

for i = 1:length(nsNeuralInfo)
    nsNeuralInfo(i).SourceEntityID = SourceEntityID(i);
    nsNeuralInfo(i).SourceUnitID = SourceUnitID(i);
end
