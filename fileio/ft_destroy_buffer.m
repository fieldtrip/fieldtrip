function ft_destroy_buffer

% FT_DESTROY_BUFFER stops the thread with the TCP server attached to
% the local Matlab instance and removes all data from memory.
%
% Use as
%   ft_destroy_buffer
%
% See also FT_CREATE_BUFFER

% clearing the mex file from memory will cause the function registered with
% mexAtExit to be executed. This function will then stop the threads and
% release all allocated memory to the operating system

clear buffer
