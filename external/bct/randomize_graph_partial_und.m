function A = randomize_graph_partial_und(A,B,maxswap)
%RANDOMIZE_GRAPH_PARTIAL_UND    swap edges/preserve degree sequence
%
%   A = RANDOMIZE_GRAPH_PARTIAL_UND(A,B,MAXSWAP) takes adjacency matrices A 
%   and B and attempts to randomize matrix A by performing MAXSWAP 
%   rewirings. The rewirings will avoid any spots where matrix B is 
%   nonzero.
%
%   Inputs:       A,      undirected adjacency matrix
%                 B,      edges to avoid
%           MAXSWAP,      number of rewirings
%
%   Outputs:      A,      randomized matrix
%
%   Richard Betzel, Indiana University, 2013
%
%   Notes:
%   1. Based on the script randmio_und.m.
%   2. Graph may become disconnected as a result of rewiring. Always
%      important to check.
%   3. A can be weighted, though the weighted degree sequence will not be
%      preserved.
%

[i,j] = find(triu(A,1));
m = length(i);
nswap = 0;
while nswap < maxswap
    while 1
        e1 = randi(m); e2 = randi(m);
        while e2 == e1
            e2 = randi(m);
        end
        a = i(e1); b = j(e1);
        c = i(e2); d = j(e2);
        if all(a~=[c,d]) && all(b~=[c,d])
            break
        end
    end
    if rand > 0.5
        i(e2) = d; j(e2) = c;
        c = i(e2); d = j(e2);
    end
    if ~(A(a,d) || A(c,b) || B(a,d) || B(c,b))
        A(a,d) = A(a,b); A(a,b) = 0;
        A(d,a) = A(b,a); A(b,a) = 0;
        A(c,b) = A(c,d); A(c,d) = 0;
        A(b,c) = A(d,c); A(d,c) = 0;
        j(e1) = d;
        j(e2) = b;
        nswap = nswap + 1;
    end
end
