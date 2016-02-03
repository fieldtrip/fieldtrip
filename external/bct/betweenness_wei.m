function BC=betweenness_wei(G)
%BETWEENNESS_WEI    Node betweenness centrality
%
%   BC = betweenness_wei(L);
%
%   Node betweenness centrality is the fraction of all shortest paths in 
%   the network that contain a given node. Nodes with high values of 
%   betweenness centrality participate in a large number of shortest paths.
%
%   Input:      L,      Directed/undirected connection-length matrix.
%
%   Output:     BC,     node betweenness centrality vector.
%
%   Notes:
%       The input matrix must be a connection-length matrix, typically
%   obtained via a mapping from weight to length. For instance, in a
%   weighted correlation network higher correlations are more naturally
%   interpreted as shorter distances and the input matrix should
%   consequently be some inverse of the connectivity matrix. 
%       Betweenness centrality may be normalised to the range [0,1] as
%   BC/[(N-1)(N-2)], where N is the number of nodes in the network.
%
%   Reference: Brandes (2001) J Math Sociol 25:163-177.
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2012

n=length(G);
% E=find(G); G(E)=1./G(E);        %invert weights
BC=zeros(n,1);                  %vertex betweenness

for u=1:n
    D=inf(1,n); D(u)=0;         %distance from u
    NP=zeros(1,n); NP(u)=1;     %number of paths from u
    S=true(1,n);                %distance permanence (true is temporary)
    P=false(n);                 %predecessors
    Q=zeros(1,n); q=n;          %order of non-increasing distance

    G1=G;
    V=u;
    while 1
        S(V)=0;                 %distance u->V is now permanent
        G1(:,V)=0;              %no in-edges as already shortest
        for v=V
            Q(q)=v; q=q-1;
            W=find(G1(v,:));                %neighbours of v
            for w=W
                Duw=D(v)+G1(v,w);           %path length to be tested
                if Duw<D(w)                 %if new u->w shorter than old
                    D(w)=Duw;
                    NP(w)=NP(v);            %NP(u->w) = NP of new path
                    P(w,:)=0;
                    P(w,v)=1;               %v is the only predecessor
                elseif Duw==D(w)            %if new u->w equal to old
                    NP(w)=NP(w)+NP(v);      %NP(u->w) sum of old and new
                    P(w,v)=1;               %v is also a predecessor
                end
            end
        end

        minD=min(D(S));
        if isempty(minD), break             %all nodes reached, or
        elseif isinf(minD),                 %...some cannot be reached:
            Q(1:q)=find(isinf(D)); break	%...these are first-in-line
        end
        V=find(D==minD);
    end

    DP=zeros(n,1);                          %dependency
    for w=Q(1:n-1)
        BC(w)=BC(w)+DP(w);
        for v=find(P(w,:))
            DP(v)=DP(v)+(1+DP(w)).*NP(v)./NP(w);
        end
    end
end