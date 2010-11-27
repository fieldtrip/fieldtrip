function ns_RESULT = ns_CloseFile(hFile);

%ns_CloseFile   Closes a neural data file
%
%   Usage:
%      ns_RESULT = ns_CloseFile(hFile)
%
%   Description:
%       Closes a previously opened file specified by the file handle hFile.
%
%   Parameters:
%       hFile	    Handle/Indentification number to an open file.
%
%   Return Values:
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_BADFILE	Invalid file handle passed to function.
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 2/4/2003

ns_RESULT = mexprog(14, hFile);