function [fid] = fiff_start_file(name)
%
% [fid] = fiff_start_file(name)
% 
% Opens a fif file for writing and writes the compulsory header tags
%
%     name           The name of the file to open. It is recommended
%                    that the name ends with .fif
%
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%

me='MNE:fiff_start_file';
if nargin ~= 1
    error(me,'File name required as an argument');
end

[fid,message] = fopen(name,'w+','ieee-be');
if fid < 0
    error(me,message);
end
%
%   Write the compulsory items
%
FIFF_FILE_ID=100;
FIFF_DIR_POINTER=101;
FIFF_FREE_LIST=106;

fiff_write_id(fid,FIFF_FILE_ID);
fiff_write_int(fid,FIFF_DIR_POINTER,-1);
fiff_write_int(fid,FIFF_FREE_LIST,-1);
%
%   Ready for more
%
return;

