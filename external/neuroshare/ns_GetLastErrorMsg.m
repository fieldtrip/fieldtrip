function [ns_RESULT, LastError] = ns_GetLastErrorMsg();

%ns_GetLastErrorMsg   Retrieves the extended error message for the last
%error
%
%   Usage:
%      [ns_RESULT, LastError] = ns_GetLastErrorMsg
%
%   Description:
%       Returns the last error message in text form to LastError. This
%       function should be called immediately following a function whose
%       return value indicates that an error occurred. The maximum size of
%       the error message text is 256 characters.
%
%   Return Values:
%       LastError	Variable to receive the extended last error message.
%       ns_RESULT   This function returns ns_OK if the file is successfully
%                   opened. Otherwise one of the following error codes is 
%                   generated:
%
%                       ns_LIBERROR	Library error
%
%   Copyright (C) 2003 Neuroshare Project
%   Author: Almut Branner
%   Last modification: 2/21/2003

[ns_RESULT, LastError] = mexprog(17);
