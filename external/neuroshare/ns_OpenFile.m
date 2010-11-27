function [ns_RESULT, hFile] = ns_OpenFile(filename);

%ns_OpenFile   Opens a neural data file
%
%   Usage:
%       [ns_RESULT, hFile] = ns_OpenFile('filename.ext')
%
%   Description:
%       Opens the file specified by filename and returns a file handle, 
%       hFile, that is used to access the opened file.
%
%   Parameters:
%       filename	Pointer to a null-terminated string that specifies the
%                   name of the file to open. 
%
%   Return Values:
%       hFile	    Handle/Indentification number to the opened file.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_TYPEERROR	Library unable to open file type
%                       ns_FILEERROR	File access or read error 
%
%   Remarks:
%       All files are opened for read-only, as no writing capabilities have
%       been implemented.  If the command succeeds in opening the file, the
%       application should call ns_CloseFile for each open file before 
%       terminating.
%
%       This function has to be called before any other Neuroshare function
%       is called. 
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 6/20/2003

[ns_RESULT, hFile] = mexprog(1, filename);