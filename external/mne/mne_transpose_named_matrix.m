function [res] = mne_transpose_named_matrix(mat)
%
% [res] = mne_transpose_named_matrix(mat)
%
% Transpose a named matrix
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/18 23:21:22  msh
%   Added mne_transform_source_space_to and mne_transpose_named_matrix
%
%

res.nrow = mat.ncol;
res.ncol = mat.nrow;
res.row_names = mat.col_names;
res.col_names = mat.row_names;
res.data      = mat.data';

return;

end

