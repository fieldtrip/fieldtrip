function fiff_write_double_complex_matrix(fid,kind,mat)
%
% fiff_write_double_complex_matrix(fid,kind,mat)
%
% Writes a double-precision complex matrix tag
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
%   Revision 1.1  2006/09/23 14:43:56  msh
%   Added routines for writing complex and double complex matrices.
%   Added routine for writing double-precision real matrix.
%
%

me='MNE:fiff_write_double_complex_matrix';

if nargin ~= 3
        error(me,'Incorrect number of arguments');
end

FIFFT_COMPLEX_DOUBLE = 21;
FIFFT_MATRIX  = bitshift(1,30);
FIFFT_MATRIX_COMPLEX_DOUBLE = bitor(FIFFT_COMPLEX_DOUBLE,FIFFT_MATRIX);
FIFFV_NEXT_SEQ=0;

datasize = 2*8*numel(mat) + 4*3;

count = fwrite(fid,int32(kind),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_MATRIX_COMPLEX_DOUBLE),'int32');
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
nrow = size(mat,1);
ncol = size(mat,2);
for j = 1:nrow
   for k = 1:ncol
      count = fwrite(fid,real(mat(j,k)),'double');
      if count ~= 1
          error(me,'write failed');
      end
      count = fwrite(fid,imag(mat(j,k)),'double');
      if count ~= 1
          error(me,'write failed');
      end
   end
end
dims(1) = size(mat,2);
dims(2) = size(mat,1);
dims(3) = 2;
count = fwrite(fid,int32(dims),'int32');
if count ~= 3
    error(me,'write failed');
end

return;

