function [CIJ] = makerandCIJ_und(N,K)
%MAKERANDCIJ_UND        Synthetic directed random network
%
%   CIJ = makerandCIJ_und(N,K);
%
%   This function generates an undirected random network
%
%   Inputs:     N,      number of vertices
%               K,      number of edges
%
%   Output:     CIJ,    undirected random connection matrix
%
%   Note: no connections are placed on the main diagonal.
%
%
% Olaf Sporns, Indiana University, 2007/2008

ind = triu(~eye(N));
i = find(ind);
rp = randperm(length(i));
irp = i(rp);

CIJ = zeros(N);
CIJ(irp(1:K)) = 1;
CIJ = CIJ+CIJ';         % symmetrize
