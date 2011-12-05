function [kden,N,K] = density_dir(CIJ)
%DENSITY        Density
%
%   kden = density_dir(CIJ);
%   [kden,N,K] = density_dir(CIJ);
%
%   Density is the fraction of present connections to possible connections.
%
%   Input:      CIJ,    directed (weighted/binary) connection matrix
%
%   Output:     kden,   density
%               N,      number of vertices
%               K,      number of edges
%
%   Notes:  Assumes CIJ is directed and has no self-connections.
%           Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);
K = nnz(CIJ);
kden = K/(N^2-N);

