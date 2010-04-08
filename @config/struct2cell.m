function c = struct2cell(s)

% STRUCT2CELL Convert structure array to cell array by first converting
% the config objection into a struct and then using standard Matlab call.

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

c = struct2cell(struct(s));

