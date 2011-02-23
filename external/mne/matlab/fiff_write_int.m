function fiff_write_int(fid,kind,data)
%
% fiff_write_int(fid,kind,data)
% 
% Writes a 32-bit integer tag to a fif file
%
%     fid           An open fif file descriptor
%     kind          Tag kind
%     data          The integers to use as data
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

me='MNE:fiff_write_int';

if nargin ~= 3
        error(me,'Incorrect number of arguments');
end

FIFFT_INT=3;
FIFFV_NEXT_SEQ=0;
nel=numel(data);
datasize=nel*4;
count = fwrite(fid,int32(kind),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_INT),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(datasize),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFV_NEXT_SEQ),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(data),'int32');
if count ~= nel
    error(me,'write failed');
end

return;
