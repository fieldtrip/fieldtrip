function fiff_write_named_matrix(fid,kind,mat)
%
% fiff_write_named_matrix(fid,kind,mat)
% 
% Writes a named single-precision floating-point matrix
%
%     fid           An open fif file descriptor
%     kind          The tag kind to use for the data
%     mat           The data matrix
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.5  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/12 10:29:03  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.3  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%   Revision 1.2  2005/12/05 20:23:21  msh
%   Added fiff_save_evoked. Improved error handling.
%
%   Revision 1.1  2005/12/05 16:01:04  msh
%   Added an initial set of fiff writing routines.
%
%

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

me='MNE:fiff_write_named_matrix';

if nargin ~= 3
        error(me,'Incorrect number of arguments');
end

fiff_start_block(fid,FIFF.FIFFB_MNE_NAMED_MATRIX);
    fiff_write_int(fid,FIFF.FIFF_MNE_NROW,mat.nrow);
    fiff_write_int(fid,FIFF.FIFF_MNE_NCOL,mat.ncol);
    if size(mat.row_names,2) > 0
       fiff_write_name_list(fid,FIFF.FIFF_MNE_ROW_NAMES,mat.row_names);    
    end
    if size(mat.col_names,2) > 0
       fiff_write_name_list(fid,FIFF.FIFF_MNE_COL_NAMES,mat.col_names);    
    end
    fiff_write_float_matrix(fid,kind,mat.data);
fiff_end_block(fid,FIFF.FIFFB_MNE_NAMED_MATRIX);

return;

