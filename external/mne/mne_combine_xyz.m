function [comb] = mne_combine_xyz(vec)
%
% function [comb] = mne_combine_xyz(vec)
%
% Compute the three Cartesian components of a vector together
%
%
% vec         - Input row or column vector [ x1 y1 z1 ... x_n y_n z_n ]
% comb        - Output vector [x1^2+y1^2+z1^2 ... x_n^2+y_n^2+z_n^2 ]
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

me = 'MNE:mne_combine_xyz';
if nargin ~= 1
    error(me,'Wrong number of arguments');
end
if size(vec,1) > size(vec,2)
    vec = vec';
end
if size(vec,1) ~= 1 || mod(size(vec,2),3) ~= 0
    error(me,'Input must be a row or a column vector with 3N components');
end

s = mne_block_diag(vec,3);
comb = full(diag(s*s'));
if size(vec,1) > size(vec,2)
    comb = comb';
end

return;

