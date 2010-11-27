function [ns_RESULT, Index] = ns_GetIndexByTime(hFile, EntityID, Time, Flag);

%ns_GetIndexByTime   Retrieves an entity index by time
%
%   Usage:
%      [ns_RESULT, Index] = ns_GetIndexByTime(hFile, EntityID, Time, Flag)
%
%   Description:
%       Searches in the file referenced by hFile for the data item
%       identified by the index EntityID.  The flag specifies whether to
%       locate the data item that starts before or after the time Time.
%       The index of the requested data item is returned in Index.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%       EntityID	Identification number of the entity in the data file.
%       Time	    Time of the data to search for
%       Flag	    Flag specifying whether the index to be retrieved
%                   belongs to the data item occurring before or after the
%                   specified time Time.  The flags are defined:
%
%               #define ns_BEFORE 	-1	// return the data entry occuring
%                                       // before and inclusive of the time
%                                       // dTime.
%               #define ns_CLOSEST	0	// return the data entry occuring
%                                       // at or closest to the time dTime 
%               #define ns_AFTER 	+1	// return the data entry occuring
%                                       // after and inclusive of the time
%                                       //dTime.
%
%   Return Values:
%       Index	    Variable to receive the entry index.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADFILE	    Invalid file handle passed to 
%                                       function
%                       ns_BADENTITY	Invalid or inappropriate entity 
%                                       identifier specified
%                       ns_FILEERROR	File access or read error
%                       ns_BADINDEX	    Unable to find an valid index given
%                                       the search parameters
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 6/20/2003

[ns_RESULT, Index] = mexprog(15, hFile, EntityID - 1, Time, Flag);

Index = Index + 1;