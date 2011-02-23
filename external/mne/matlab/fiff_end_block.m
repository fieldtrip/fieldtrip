function fiff_end_block(fid, kind)
%
% fiff_end_block(fid, kind)
%
% Writes a FIFF_BLOCK_END tag
%
%     fid           An open fif file descriptor
%     kind          The block kind to end
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

me = 'MNE:fiff_end_block';

if nargin ~= 2
    error(me, 'Incorrect number of arguments');
end

FIFF_BLOCK_END = 105;
fiff_write_int(fid, FIFF_BLOCK_END, kind);

return
