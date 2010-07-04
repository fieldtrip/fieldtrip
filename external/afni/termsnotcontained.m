function m = termsnotcontained(t)
%TERMSNOTCONTAINED Creates a logical matrix indicating term containment
%   m(i,j)==1  iff  t(i) is not contained by t(j)
%   set diagonals to 1 because we want proper containment

m = (t*~t') > 0;
m(1:(size(t,1)+1):end) = 1;
