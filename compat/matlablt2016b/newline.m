function C = newline

% NEWLINE Create newline character 
%   C = newline creates the character representing a newline.
%   newline is equivalent to char(10) or sprintf('\n').
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

C = char(10);
