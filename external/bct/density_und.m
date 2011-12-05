function [kden,N,K] = density_und(CIJ)
%DENSITY        Density
%
%   kden = density_und(CIJ);
%   [kden,N,K] = density_und(CIJ);
%
%   Density is the fraction of present connections to possible connections.
%
%   Input:      CIJ,    undirected (weighted/binary) connection matrix
%
%   Output:     kden,   density
%               N,      number of vertices
%               K,      number of edges
%
%   Notes:  Assumes CIJ is undirected and has no self-connections.
%           Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008


% Modification history:
% 2009-10: K fixed to sum over one half of CIJ [Tony Herdman, SFU]

N = size(CIJ,1);
K = nnz(triu(CIJ));
kden = K/((N^2-N)/2);

