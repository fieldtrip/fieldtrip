function [ names ] = fiff_split_name_list(list);
%
% [names] = fiff_split_name_list(list)
%
%
% Split a name list containing colon-separated entries into a cell array
% containing the strings
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2006/04/27 22:38:37  msh
%   Splitting an empty list now results in an empty output.
%   Added fiff_write_double and fiff_write_short
%   Write an array of floats, ints, and doubles instead of just one value.
%   Fixed help text of fiff_pick_channels.
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

rem = list;
nnames = 0;
if isempty(rem)
   names = [];
   return;
end
while true
    [ str, rem ]  = strtok(rem,':');
    if isempty(str)
        break;
    end
    nnames = nnames + 1;
    names{nnames} = str;
end

return;

