function result = is_orthogonal_matrix(P)
%
% Orthogonal Matrices
% 
% is_orthogonal_matrix(P) determines if the matrix P is an orthogonal
% matrix. An error is returned if a matrix that is not square is attempted
% to be determined for orthogonality.
%
%  Function written by Anthony Russo, downloaded from MatlabCentral.
%

matrix_size = size(P);

m = matrix_size(1,1);
n = matrix_size(1,2);

tolerance = 10^-10;

if m ~= n
    error('Only square matrices can be orthogonal.');
else
    count = 0;

    identity_matrix = P*P';

    if det(P) ~= 0
        for i = 1:m
            if abs(identity_matrix(i,i) - 1) <= tolerance
                count = count + 1;
            else
                break
            end
        end
    end

    if count == m
        result = 1;
    else
        result = 0;
    end
end
