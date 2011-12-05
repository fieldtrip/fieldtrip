function D=distance_wei(G)
%DISTANCE_WEI       Distance matrix
%
%   D = distance_wei(W);
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path 
%   from node u to node v. The average shortest path length is the 
%   characteristic path length of the network.
%
%   Input:      W,      weighted directed/undirected connection matrix
%
%   Output:     D,      distance matrix
%
%   Notes:
%       The input matrix must be a mapping from weight to distance. For 
%   instance, in a weighted correlation network, higher correlations are 
%   more naturally interpreted as shorter distances, and the input matrix 
%   should consequently be some inverse of the connectivity matrix.
%       Lengths between disconnected nodes are set to Inf.
%       Lengths on the main diagonal are set to 0.
%
%   Algorithm: Dijkstra's algorithm.
%
%
%   Mika Rubinov, UNSW, 2007-2010.

%Modification history
%2007: original
%2009-08-04: min() function vectorized

n=length(G);
D=zeros(n); D(~eye(n))=inf;                 %distance matrix

for u=1:n
    S=true(1,n);                            %distance permanence (true is temporary)
    G1=G;
    V=u;
    while 1
        S(V)=0;                             %distance u->V is now permanent
        G1(:,V)=0;                          %no in-edges as already shortest
        for v=V
            W=find(G1(v,:));                %neighbours of shortest nodes
            D(u,W)=min([D(u,W);D(u,v)+G1(v,W)]); %smallest of old/new path lengths
        end

        minD=min(D(u,S));
        if isempty(minD)||isinf(minD),      %isempty: all nodes reached;
            break,                          %isinf: some nodes cannot be reached
        end;

        V=find(D(u,:)==minD);
    end
end