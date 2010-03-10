function destroy_buffer

% DESTROY_BUFFER stops the thread with the TCP server attached to
% the local Matlab instance and removes all data from memory.
%
% Use as
%   destroy_buffer
%
% See also CREATE_BUFFER

% clearing the mex file from memory will cause the function registered with
% mexAtExit to be executed. This function will then stop the threads and
% release all allocated memory to the operating system

clear buffer
