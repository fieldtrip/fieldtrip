function [CIJ] = makerandCIJ_dir(N,K)
%MAKERANDCIJ_DIR        Synthetic directed random network
%
%   CIJ = makerandCIJ_dir(N,K);
%
%   This function generates a directed random network
%
%   Inputs:     N,      number of vertices
%               K,      number of edges
%
%   Output:     CIJ,    directed random connection matrix
%
%   Note: no connections are placed on the main diagonal.
%
%
% Olaf Sporns, Indiana University, 2007/2008

ind = ~eye(N);
i = find(ind);
rp = randperm(length(i));
irp = i(rp);

CIJ = zeros(N);
CIJ(irp(1:K)) = 1;
