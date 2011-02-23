function fiff_end_file(fid)
%
% fiff_end_file(fid)
%
% Writes the closing tags to a fif file and closes the file
%
%     fid           An open fif file descriptor
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

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me = 'MNE:fiff_end_file';

if nargin ~= 1
    error(me, 'An open file id required as an argument');
end

datasize = 0;
count = fwrite(fid, int32(FIFF.FIFF_NOP), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(FIFF.FIFFT_VOID), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(datasize), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(FIFF.FIFFV_NEXT_NONE), 'int32');
if count ~= 1
    error(me, 'write failed');
end
fclose(fid);

return

