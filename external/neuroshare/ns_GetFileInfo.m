function [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);

%ns_GetFileInfo   Retrieves file information and entity counts
%
%   Usage:
%       [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile)
%
%   Description:
%       Provides general information about the data file referenced by 
%       hFile. This information is returned in the structure nsFileInfo.
%
%   Parameters:
%       hFile	        Handle/Indentification number to an open file.
%
%   Return Values:
%       nsFileInfo  ns_FILEINFO structure that receives the file
%                   information.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_FILEERROR	File access or read error
%                       ns_BADFILE	    Invalid file handle passed to
%                                       function
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 2/4/2003

[ns_RESULT, nsFileInfo] = mexprog(3, hFile);
