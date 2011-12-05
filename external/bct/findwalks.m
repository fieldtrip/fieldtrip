function [Wq,twalk,wlq] = findwalks(CIJ)
%FINDWALKS      Network walks
%
%   [Wq,twalk,wlq] = findwalks(CIJ);
%
%   Walks are sequences of linked nodes, that may visit a single node more
%   than once. This function finds the number of walks of a given length, 
%   between any two nodes.
%
%   Input:      CIJ         binary (directed/undirected) connection matrix
%
%   Outputs:    Wq          3D matrix, Wq(i,j,q) is the number of walks
%                           from 'i' to 'j' of length 'q'.
%               twalk       total number of walks found
%               wlq         walk length distribution as function of 'q'
%
%   Notes: Wq grows very quickly for larger N,K,q. Weights are discarded.
%
%   Algorithm: algebraic path count
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

% ensure CIJ is binary...
CIJ = double(CIJ~=0);

N = size(CIJ,1);
Wq = zeros(N,N,N);
CIJpwr = CIJ;
Wq(:,:,1) = CIJ;
for q=2:N
   CIJpwr = CIJpwr*CIJ;
   Wq(:,:,q) = CIJpwr;
end;

% total number of walks
twalk = sum(sum(sum(Wq)));

% walk length distribution
wlq = reshape(sum(sum(Wq)),1,N);

