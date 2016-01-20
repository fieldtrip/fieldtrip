function [CIJ] = makelatticeCIJ(N,K)
%MAKELATTICECIJ     Synthetic lattice network
%
%   CIJ = makelatticeCIJ(N,K);
%
%   This function generates a directed lattice network without toroidal 
%   boundary counditions (i.e. no ring-like "wrapping around").
%
%   Inputs:     N,      number of vertices
%               K,      number of edges
%
%   Outputs:    CIJ,    connection matrix
%
%   Note: The lattice is made by placing connections as close as possible 
%   to the main diagonal, without wrapping around. No connections are made 
%   on the main diagonal. In/Outdegree is kept approx. constant at K/N.
%
%
%   Olaf Sporns, Indiana University, 2005/2007

% initialize
CIJ = zeros(N);
CIJ1 = ones(N);
KK = 0;
cnt = 0;
seq = 1:N-1;

% fill in
while (KK<K)
    cnt = cnt + 1;
    dCIJ = triu(CIJ1,seq(cnt))-triu(CIJ1,seq(cnt)+1);
    dCIJ = dCIJ+dCIJ';
    CIJ = CIJ + dCIJ;
    KK = sum(sum(CIJ));
end;

% remove excess connections
overby = KK-K;
if(overby>0)
    [i,j] = find(dCIJ);
    rp = randperm(length(i));
    for ii=1:overby
        CIJ(i(rp(ii)),j(rp(ii))) = 0;
    end;
end;
