function [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ns_GetSegmentData(hFile, EntityID, Index);

%ns_GetSegmentData   Retrieves segment data by index
%
%   Usage:
%      [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = 
%                               ns_GetSegmentData(hFile, EntityID, Index)
%
%   Description:
%       Returns the Segment data values in entry Index of the entity
%       EntityID from the file referenced by hFile. The data values are
%       returned in Data. The timestamp of the entry id returned in
%       TimeStamp. The number of samples written to Data is returned in
%       SampleCount. 
%       The data variable should be accessed as a 2-dimensional array for
%       samples and sources.  
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID    Identification number of the entity in the data file.
%       Index	    Index number of the requested Segment data item.
%
%   Remarks:
%       A zero unit ID is unclassified, then follow unit 1, 2, 3, etc. Unit
%       255 is noise.
%
%   Return Values:
%       TimeStamp	Time stamp of the requested Segment data item.
%       Data	    Variable to receive the requested data.
%       SampleCount	Number of samples returned in the data variable.
%       UnitID	    Unit classification code for the Segment Entity.
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
%   Last modification: 10/24/2003


[ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = mexprog(11, hFile, EntityID - 1, Index - 1);
Data = squeeze(Data);
if (size(Data, 2) == 1)
    Data = Data';
end;

ind = find(UnitID == 1);
UnitID(ind) = 255;

ind = find((UnitID > 1) & (UnitID < 255));
UnitID(ind) = log2(UnitID(ind));