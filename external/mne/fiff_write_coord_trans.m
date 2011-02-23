function fiff_write_coord_trans(fid,trans)
%
% fiff_write_coord_trans(fid,trans)
% 
% Writes a coordinate transformation structure
%
%     fid           An open fif file descriptor
%     trans         The coordinate transfomation structure
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.2  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%

me='MNE:fiff_write_coord_trans';

if nargin ~= 2
        error(me,'Incorrect number of arguments');
end

FIFF_COORD_TRANS=222;
FIFFT_COORD_TRANS_STRUCT=35;
FIFFV_NEXT_SEQ=0;


%?typedef struct _fiffCoordTransRec {
%  fiff_int_t   from;                   /*!< Source coordinate system. */
%  fiff_int_t   to;                     /*!< Destination coordinate system. */
%  fiff_float_t rot[3][3];              /*!< The forward transform (rotation part) */
%  fiff_float_t move[3];                /*!< The forward transform (translation part) */
%  fiff_float_t invrot[3][3];           /*!< The inverse transform (rotation part) */
%  fiff_float_t invmove[3];             /*!< The inverse transform (translation part) */
%} *fiffCoordTrans, fiffCoordTransRec;  /*!< Coordinate transformation descriptor */

datasize=4*2*12 + 4*2;
count = fwrite(fid,int32(FIFF_COORD_TRANS),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(FIFFT_COORD_TRANS_STRUCT),'int32');
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
%   Start writing fiffCoordTransRec
%
count = fwrite(fid,int32(trans.from),'int32');
if count ~= 1
    error(me,'write failed');
end
count = fwrite(fid,int32(trans.to),'int32');
if count ~= 1
    error(me,'write failed');
end
%
%   The transform...
%
rot=trans.trans(1:3,1:3)';
move=trans.trans(1:3,4)';
count = fwrite(fid,single(rot),'single');
if count ~= 9
    error(me,'write failed');
end
count = fwrite(fid,single(move),'single');
if count ~= 3
    error(me,'write failed');
end
%
%   ...and its inverse
%
trans_inv=inv(trans.trans);
rot=trans_inv(1:3,1:3)';
move=trans_inv(1:3,4)';
count = fwrite(fid,single(rot),'single');
if count ~= 9
    error(me,'write failed');
end
count = fwrite(fid,single(move),'single');
if count ~= 3
    error(me,'write failed');
end
