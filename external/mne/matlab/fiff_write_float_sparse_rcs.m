function fiff_write_float_sparse_rcs(fid,kind,mat)
%
% fiff_write_float_sparsce_rcs(fid,kind,mat)
% 
% Writes a single-precision sparse (RCS) floating-point matrix tag
%
%     fid           An open fif file descriptor
%     kind          The tag kind
%     mat           The data matrix
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%

me='MNE:fiff_write_float_sparse_rcs';

if nargin ~= 3
   error(me,'Incorrect number of arguments');
end

if ~issparse(mat)
   error(me,'Input should be a sparse matrix');
end

if length(size(mat)) ~= 2
   error(me,'Input should be a two-dimensional matrix');
end

FIFFT_FLOAT  = 4;
FIFFT_MATRIX = bitshift(16416,16);    % 4020
FIFFT_MATRIX_FLOAT_RCS = bitor(FIFFT_FLOAT,FIFFT_MATRIX);
FIFFV_NEXT_SEQ=0;
%
%    nnz values
%    nnz column indices
%    nrow+1 pointers
%    dims
%    nnz
%    ndim
%
nnzm = nnz(mat);
nrow = size(mat,1);
datasize = 4*nnzm + 4*nnzm + 4*(nrow+1) + 4*4;
%
%    Nonzero entries
%
[ s(:,1), s(:,2), s(:,3) ] = find(mat);
s = sortrows(s,1);
[ rows, starts ] = unique(s(:,1),'first');

count = fwrite(fid,int32(kind),'int32');
if count ~= 1
   error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_MATRIX_FLOAT_RCS),'int32');
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
%
%   The data values
%
count = fwrite(fid,single(s(:,3)),'single');
if count ~= nnzm
   error(me,'write failed');
end
%
%   Column indices
%
count = fwrite(fid,int32(s(:,2)-1),'int32');
if count ~= nnzm
   error(me,'write failed');
end
%
%   Pointers
%
ptrs = -ones(1,nrow+1);
for k = 1:length(rows)
   ptrs(rows(k)) = starts(k) - 1;
end
ptrs(nrow+1) = nnzm;
%
%   Fill in pointers for empty rows
%
for k = nrow:-1:1
   if ptrs(k) < 0
      ptrs(k) = ptrs(k+1);
   end
end
%
count = fwrite(fid,int32(ptrs),'int32');
if count ~= nrow+1
   error(me,'write failed');
end
%
%   Dimensions
%
dims(1) = nnz(mat);
dims(2) = size(mat,1);
dims(3) = size(mat,2);
dims(4) = 2;
count = fwrite(fid,int32(dims),'int32');
if count ~= 4
   error(me,'write failed');
end

return;

