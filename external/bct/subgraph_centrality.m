function Cs = subgraph_centrality(CIJ)
% SUBGRAPH_CENTRALITY     subgraph centrality of a network
%
%   Cs = subgraph_centrality(CIJ)
%
%   The subgraph centrality of a node is a weighted sum of closed walks of
%   different lengths in the network starting and ending at the node. This
%   function returns a vector of subgraph centralities for each node of the
%   network.
%
%   Inputs:     CIJ,        adjacency matrix (binary)
%
%   Outputs:     Cs,        subgraph centrality
%
%   Reference: Estrada and Rodriguez-Velasquez (2005) Phys Rev E 71, 056103
%              Estrada and Higham (2010) SIAM Rev 52, 696.
%
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012

[V,lambda] = eig(CIJ);                 % Compute the eigenvectors and
lambda     = diag(lambda);             % eigenvalues.
V2         = V.^2;                     % Matrix of squares of the eigenvectors elements.
Cs         = real(V2 * exp(lambda));   % Compute eigenvector centrality. Lop off imaginary part remaining due to precision error.