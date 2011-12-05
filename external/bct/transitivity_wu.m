function T=transitivity_wu(W)
%TRANSITIVITY_WU    Transitivity
%
%   T = transitivity_wu(W);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      W       weighted undirected connection matrix
%
%   Output:     T       transitivity scalar
%
%   Reference: Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW, 2010

K=sum(W~=0,2);            	
cyc3=diag((W.^(1/3))^3);           
T=sum(cyc3)./sum((K.*(K-1)));       %transitivity
