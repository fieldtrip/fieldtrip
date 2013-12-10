function E=efficiency(G,local)
%EFFICIENCY     Global efficiency, local efficiency.
%
%   Eglob = efficiency(A);
%   Eloc = efficiency(A,1);
%
%   The global efficiency is the average of inverse shortest path length, 
%   and is inversely related to the characteristic path length.
%
%   The local efficiency is the global efficiency computed on the
%   neighborhood of the node, and is related to the clustering coefficient.
%
%   Inputs:     A,              binary undirected connection matrix
%               local,          optional argument
%                               (local=1 computes local efficiency)
%
%   Output:     Eglob,          global efficiency (scalar)
%               Eloc,           local efficiency (vector)
%
%
%   Algorithm: algebraic path count
%
%   Reference: Latora and Marchiori (2001) Phys Rev Lett 87:198701.
%
%
%   Mika Rubinov, UNSW, 2008-2010

if ~exist('local','var')
    local=0;
end

if local                                %local efficiency
    N=length(G);                        %number of nodes
    E=zeros(N,1);                       %local efficiency

    for u=1:N
        V=find(G(u,:));                 %neighbors
        k=length(V);                    %degree
        if k>=2;                        %degree must be at least two
            e=distance_inv(G(V,V));
            E(u)=sum(e(:))./(k^2-k);	%local efficiency
        end
    end
else
    N=length(G);
    e=distance_inv(G);
    E=sum(e(:))./(N^2-N);               %global efficiency
end


function D=distance_inv(g)
D=eye(length(g));
n=1;
nPATH=g;                        %n-path matrix
L=(nPATH~=0);                   %shortest n-path matrix

while find(L,1);
    D=D+n.*L;
    n=n+1;
    nPATH=nPATH*g;
    L=(nPATH~=0).*(D==0);
end

D(~D)=inf;                      %disconnected nodes are assigned d=inf;
D=1./D;                         %invert distance
D=D-eye(length(g));