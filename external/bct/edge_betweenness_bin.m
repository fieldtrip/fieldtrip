function [EBC,BC]=edge_betweenness_bin(G)
%EDGE_BETWEENNESS_BIN    Edge betweenness centrality
%
%   EBC = edge_betweenness_bin(A);
%   [EBC BC] = edge_betweenness_bin(A);
%
%   Edge betweenness centrality is the fraction of all shortest paths in 
%   the network that contain a given edge. Edges with high values of 
%   betweenness centrality participate in a large number of shortest paths.
%
%   Input:      A,      binary (directed/undirected) connection matrix.
%
%   Output:     EBC,    edge betweenness centrality matrix.
%               BC,     node betweenness centrality vector.
%
%   Note: Betweenness centrality may be normalised to the range [0,1] as
%   BC/[(N-1)(N-2)], where N is the number of nodes in the network.
%
%   Reference: Brandes (2001) J Math Sociol 25:163-177.
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2012


n=length(G);
BC=zeros(n,1);                  %vertex betweenness
EBC=zeros(n);                   %edge betweenness

for u=1:n
    D=false(1,n); D(u)=1;      	%distance from u
    NP=zeros(1,n); NP(u)=1;     %number of paths from u
    P=false(n);                 %predecessors
    Q=zeros(1,n); q=n;          %order of non-increasing distance

    Gu=G;
    V=u;
    while V
        Gu(:,V)=0;              %remove remaining in-edges
        for v=V
            Q(q)=v; q=q-1;
            W=find(Gu(v,:));                %neighbours of v
            for w=W
                if D(w)
                    NP(w)=NP(w)+NP(v);      %NP(u->w) sum of old and new
                    P(w,v)=1;               %v is a predecessor
                else
                    D(w)=1;
                    NP(w)=NP(v);            %NP(u->w) = NP of new path
                    P(w,v)=1;               %v is a predecessor
                end
            end
        end
        V=find(any(Gu(V,:),1));
    end
    if ~all(D)                              %if some vertices unreachable,
        Q(1:q)=find(~D);                    %...these are first-in-line
    end

    DP=zeros(n,1);                          %dependency
    for w=Q(1:n-1)
        BC(w)=BC(w)+DP(w);
        for v=find(P(w,:))
            DPvw=(1+DP(w)).*NP(v)./NP(w);
            DP(v)=DP(v)+DPvw;
            EBC(v,w)=EBC(v,w)+DPvw;
        end
    end
end