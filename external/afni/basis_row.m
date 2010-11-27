function basis = basis_row(A)
%
% Bases
% 
% basis_row(A) produces a basis for the subspace of Eucldiean n-space 
% spanned by the vectors {u1,u2,...}, where the matrix A is formed from 
% these vectors as its rows. That is, the subspace is the row space of A. 
% The rows of the matrix that is returned are the basis vectors for the 
% subspace. In general, these basis vectors will not be a subset of the 
% original vectors. An error is returned if a basis for the zero vector 
% space is attempted to be produced.
%
% For example, if the vector space V = span{u1,u2,...}, where u1,u2,... are
% row vectors, then set A to be [u1;u2;...].
%
% For example, if the vector space V = Row(B), where B is an m x n matrix,
% then set A to be equal to B.
%
%  Function written by Anthony Russo, downloaded from MatlabCentral.
%

matrix_size = size(A);

m = matrix_size(1,1);
n = matrix_size(1,2);

if A == zeros(m,n)
    error('There does not exist a basis for the zero vector space.');
elseif m == 1
    basis = A;
else
    flag = 0;

    if m == 2
        multiple = A(2,1)/A(1,1);
        count = 0;

        for i = 1:n
            if A(2,i)/A(1,i) == multiple
                count = count + 1;
            end
        end

        if count == n
            basis = A(1,1:n);
            flag = 1;
        end
    end

    if flag == 0
        A = rref(A);

        for i = 1:m
            if A(i,1:n) == zeros(1,n)
                break
            else 
                basis(i,1:n) = A(i,1:n);
            end
        end
    end
end
