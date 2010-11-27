function [ns_RESULT, ContCount, Data] = ns_GetAnalogData(hFile, EntityID, StartIndex, IndexCount);

%ns_GetAnalogData   Retrieves analog data by index
%
%   Usage:
%      [ns_RESULT, ContCount, Data] =  
%               ns_GetAnalogData(hFile, EntityID, StartIndex, IndexCount)
%   
%   Description:
%       Returns the data values associated with the Analog Entity indexed
%       EntityID in the file referenced by hFile.  The index of the first
%       data value is StartIndex and the requested number of data samples
%       is given by IndexCount.  The requested data values are returned
%       in the variable Data. 
%       Although the samples in an analog entity are indexed, they are not
%       guaranteed to be continuous in time and may contain gaps between
%       some of the indexes.  When the requested data is returned,
%       ContCount contains the number of continuous data points present
%       in the data (starting at StartIndex).
%       If the index range specified by StartIndex and IndexCount contains
%       invalid indexes, the function will return ns_BADINDEX.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%       EntityID	    Identification number of the Analog Entity in the
%                       data file.
%       StartIndex	    Starting index number of the analog data item.
%       IndexCount	    Number of analog values to retrieve.
%
%   Return Values:
%       ContCount	    Number of continuous data values starting with
%                       StartIndex.  This field is ignored if the pointer
%                       is set to NULL.
%       Data	        Array of double precision values to receive the
%                       analog data.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADENTITY	Invalid or inappropriate entity 
%                                       identifier specified
%                       ns_BADINDEX	    Invalid entity index or range 
%                                       specified
%                       ns_FILEERROR	File access or read error
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 6/20/2003

[ns_RESULT, ContCount, Data] = mexprog(8, hFile, EntityID - 1, StartIndex - 1, IndexCount);
