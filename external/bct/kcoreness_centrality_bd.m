function  [coreness,kn] = kcoreness_centrality_bd(CIJ)
%KCORENESS_CENTRALITY_BD       K-coreness centrality
%
%   [coreness,kn] = kcoreness_centrality_bd(CIJ)
%
%   The k-core is the largest subgraph comprising nodes of degree at least
%   k. The coreness of a node is k if the node belongs to the k-core but
%   not to the (k+1)-core. This function computes k-coreness of all nodes
%   for a given binary directed connection matrix.
%
%   input:          CIJ,        connection/adjacency matrix (binary, directed)
%
%   output:    coreness,        node coreness.
%                    kn,        size of k-core
%
%   References: e.g. Hagmann et al. (2008) PLoS Biology
%
%   Olaf Sporns, Indiana University, 2007/2008/2010/2012

N = size(CIJ,1);

coreness = zeros(1,N); kn = zeros(1,N);
for k=1:N
    [CIJkcore,kn(k)] = kcore_bd(CIJ,k);
    ss = sum(CIJkcore)>0;
    coreness(ss) = k;
end;
