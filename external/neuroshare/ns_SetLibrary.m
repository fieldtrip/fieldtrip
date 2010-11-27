function [ns_RESULT] = ns_SetLibrary(filename);

%ns_SetDLL   Opens a Neuroshare Shared Library (.DLL or .so) in prepearation for other work
%
%   Usage:
%       [ns_RESULT] = ns_SetLibrary('filename.dll')
%
%   Description:
%       Opens the dynamic linked library specified by filename
%
%   Parameters:
%       filename	Pointer to a null-terminated string that specifies the
%                   name of the file to open. 
%
%   Return Values:
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_TYPEERROR	Library unable to open file type
%                       ns_LIBERROR	    File access or read error 
%
%   Remarks:
%       Any already opened libraries, as well as any open files will be closed
%       immediately.
%
%       All files are opened for read-only, as no writing capabilities have
%       been implemented.  If the command succeeds in opening the file, the
%       application should call ns_CloseFile for each open file before 
%       terminating.
%
%       This function has to be called before any other Neuroshare function
%       is called. 
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Kirk Korver
%   Last modification: 12/17/2003

[ns_RESULT] = mexprog(18, filename);