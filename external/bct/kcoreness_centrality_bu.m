function  [coreness,kn] = kcoreness_centrality_bu(CIJ)
%KCORENESS_CENTRALITY_BU       K-coreness centrality
%
%   [coreness,kn] = kcoreness_centrality_bu(CIJ)
%
%   The k-core is the largest subgraph comprising nodes of degree at least
%   k. The coreness of a node is k if the node belongs to the k-core but
%   not to the (k+1)-core. This function computes the coreness of all nodes
%   for a given binary undirected connection matrix.
%
%   input:          CIJ,        connection/adjacency matrix (binary, undirected)
%
%   output:    coreness,        node coreness.
%                    kn,        size of k-core
%
%   References: e.g. Hagmann et al. (2008) PLoS Biology
%
%   Olaf Sporns, Indiana University, 2007/2008/2010/2012

N = size(CIJ,1);

% determine if the network is undirected - if not, compute coreness on the
% corresponding undirected network
CIJund = CIJ+CIJ';
if (any(CIJund(:)>1))
    CIJ = double(CIJund>0);
end;

coreness = zeros(1,N); kn = zeros(1,N);
for k=1:N
    [CIJkcore,kn(k)] = kcore_bu(CIJ,k);
    ss = sum(CIJkcore)>0;
    coreness(ss) = k;
end;
