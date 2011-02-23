function fiff_write_dau16(fid, kind, data)
%
% fiff_write_dau16(fid, kind, data)
%
% Writes a 16-bit integer tag to a fif file
%
%     fid           An open fif file descriptor
%     kind          Tag kind
%     data          The integers to use as data
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   No part of this program may be photocopied, reproduced,
%   or translated to another program language without the
%   prior written consent of the author.
%

me = 'MNE:fiff_write_dau16';

if nargin ~= 3
    error(me, 'Incorrect number of arguments');
end

FIFFT_DAU_PACK16 = 16;
FIFFV_NEXT_SEQ = 0;
nel = numel(data);
datasize = nel * 2;
count = fwrite(fid, int32(kind), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(FIFFT_DAU_PACK16), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(datasize), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int32(FIFFV_NEXT_SEQ), 'int32');
if count ~= 1
    error(me, 'write failed');
end
count = fwrite(fid, int16(data), 'int16');
if count ~= nel
    error(me, 'write failed');
end

return
