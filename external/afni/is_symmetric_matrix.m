function result = is_symmetric_matrix(A)
%
% Symmetric Matrices
% 
% is_symmetric_matrix(A) determines if the matrix A is a symmetric
% matrix. An error is returned if a matrix that is not square is attempted
% to be determined for symmetry.
%
%  Function written by Anthony Russo, downloaded from MatlabCentral.
%

matrix_size = size(A);

m = matrix_size(1,1);
n = matrix_size(1,2);

if m ~= n
    error('Only square matrices can be symmetric.');
else
    if A == A'
        result = 1;
    else
        result = 0;
    end
end
