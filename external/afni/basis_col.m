function basis = basis_col(A)
%
% Bases
% 
% basis_col(A) produces a basis for the subspace of Eucldiean n-space 
% spanned by the vectors {u1,u2,...}, where the matrix A is formed from 
% these vectors as its columns. That is, the subspace is the column space 
% of A. The columns of the matrix that is returned are the basis vectors 
% for the subspace. These basis vectors will be a subset of the original 
% vectors. An error is returned if a basis for the zero vector space is 
% attempted to be produced.
%
% For example, if the vector space V = span{u1,u2,...}, where u1,u2,... are
% row vectors, then set A to be [u1' u2' ...].
%
% For example, if the vector space V = Col(B), where B is an m x n matrix,
% then set A to be equal to B.
%
%  Function written by Anthony Russo, downloaded from MatlabCentral.
%

matrix_size = size(A);

m = matrix_size(1,1);
n = matrix_size(1,2);

if A == zeros(m,n)
    error('There does not exist a basis for the zero vector space.');
elseif n == 1
    basis = A;
else
    flag = 0;

    if n == 2
        multiple = A(1,2)/A(1,1);
        count = 0;
 
        for i = 1:m
            if A(i,2)/A(i,1) == multiple
                count = count + 1;
            end
        end
 
        if count == m
            basis = A(1:m,1);
            flag = 1;
        end
    end

    if flag == 0
        [ref_A pivot_columns] = ref(A);

        for i = 1:size(pivot_columns,2)
            B(1:m,i) = A(1:m,pivot_columns(1,i));
        end

        basis = B;
    end
end
