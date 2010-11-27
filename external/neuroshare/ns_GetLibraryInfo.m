function [ns_RESULT, nsLibraryInfo] = ns_GetLibraryInfo();

%ns_GetLibraryInfo   Get library version information
%
%   Usage:
%       [ns_RESULT, nsLibraryInfo] = ns_GetLibraryInfo()
%
%   Description:
%       Obtains information about the API library.
%
%   Return Values:
%       nsLibraryInfo   Structure to receive library version information.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_LIBERROR	Library Error
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 8/11/2003

[ns_RESULT, nsLibraryInfo] = mexprog(2);
