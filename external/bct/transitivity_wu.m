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
%   Note:      All weights must be between 0 and 1.
%              This may be achieved using the weight_conversion.m function,
%              W_nrm = weight_conversion(W, 'normalize');
%
%   Reference: Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
%              based on Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2010-2015

%   Modification history:
%   2010: Original
%   2015: Expanded documentation

K    = sum(W~=0,2);            	
cyc3 = diag((W.^(1/3))^3);           
T    = sum(cyc3)./sum((K.*(K-1)));       %transitivity
