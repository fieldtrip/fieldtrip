function unit_vector = normalise(v)
%
% Normalisation
% 
% normalise(v) produces a unit vector that is in the same direction as the 
% vector v. Both column vectors and row vectors can be normalised. An error
% is returned if the zero vector is attempted to be normalised.
%
% For example, normalise([1 1]) produces the vector [1/sqrt(2) 1/sqrt(2)].
% Notice that this new vector is of length 1 whilst remaining in the same
% direction as the vector [1 1].
%
%  Function written by Anthony Russo, downloaded from MatlabCentral.
%

vector_size = size(v);

m = vector_size(1,1);
n = vector_size(1,2);

if v == zeros(m,n)
    error('The zero vector cannot be normalised.');
else
    unit_vector = v/norm(v);
end
