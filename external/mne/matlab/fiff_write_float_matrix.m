function fiff_write_float_matrix(fid,kind,mat)
%
% fiff_write_float_matrix(fid,kind,mat)
% 
% Writes a single-precision floating-point matrix tag
%
%     fid           An open fif file descriptor
%     kind          The tag kind
%     mat           The data matrix
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

me='MNE:fiff_write_float_matrix';

if nargin ~= 3
   error(me,'Incorrect number of arguments');
end

if length(size(mat)) ~= 2
   error(me,'Input should be a two-dimensional matrix');
end

FIFFT_FLOAT  = 4;
FIFFT_MATRIX = bitshift(1,30);
FIFFT_MATRIX_FLOAT = bitor(FIFFT_FLOAT,FIFFT_MATRIX);
FIFFV_NEXT_SEQ=0;

datasize = 4*numel(mat) + 4*3;

count = fwrite(fid,int32(kind),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_MATRIX_FLOAT),'int32');
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
count = fwrite(fid,single(mat'),'single');
if count ~= numel(mat)
    error(me,'write failed');
end
dims(1) = size(mat,2);
dims(2) = size(mat,1);
dims(3) = 2;
count = fwrite(fid,int32(dims),'int32');
if count ~= 3
    error(me,'write failed');
end

return;

