function [C_tri]=transitivity_bu(A)
%TRANSITIVITY_BU    Transitivity
%
%   T = transitivity_bu(A);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      A       binary undirected connection matrix
%
%   Output:     T       transitivity scalar
%
%   Reference: e.g. Humphries et al. (2008) Plos ONE 3: e0002051
%
%
%   Alexandros Goulas, Maastricht University, 2010

    C_tri = trace(A^3) / (sum(sum(A^2)) - trace(A^2));

return;