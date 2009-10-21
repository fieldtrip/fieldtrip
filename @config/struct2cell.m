function c = struct2cell(s)

% STRUCT2CELL Convert structure array to cell array by first converting
% the config objection into a struct and then using standard Matlab call.

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: struct2cell.m,v $
% Revision 1.1  2008/12/09 08:49:02  roboos
% new function
%

c = struct2cell(struct(s));

