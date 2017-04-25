function [IA, IB] = match_val(A, B)

A = A(:);
B = B(:);

[C,IA,IB] = intersect(A,B,'rows','stable');