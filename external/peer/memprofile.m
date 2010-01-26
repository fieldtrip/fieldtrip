function varargout = memprofile(varargin)

%  MEMPROFILE ON starts the profiler and clears previously recorded
%  profile statistics.
% 
%  MEMPROFILE OFF stops the profiler.
% 
%  MEMPROFILE RESUME restarts the profiler without clearing
%  previously recorded memory statistics.
% 
%  MEMPROFILE CLEAR clears all recorded profile statistics.
% 
%  STATS = MEMPROFILE('INFO') suspends the profiler and returns
%  a structure containing the current profiler statistics.

error('cannot locate mex file');

