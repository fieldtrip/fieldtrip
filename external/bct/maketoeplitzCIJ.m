function  [CIJ] = maketoeplitzCIJ(N,K,s)
%MAKETOEPLITZCIJ    A synthetic directed network with Gaussian drop-off of
%                   connectivity with distance
%
%   CIJ = maketoeprandCIJ(N,K,s)
%
%   This function generates a directed network with a Gaussian drop-off in
%   edge density with increasing distance from the main diagonal. There are
%   toroidal boundary counditions (i.e. no ring-like "wrapping around").
%
%   Inputs:     N,      number of vertices
%               K,      number of edges
%               s,      standard deviation of toeplitz
%
%   Output:     CIJ,    connection matrix
%
%   Note: no connections are placed on the main diagonal.
%
%
% Olaf Sporns, Indiana University, 2005/2007

profile = normpdf(1:N-1,0.5,s);
template = toeplitz([0 profile],[0 profile]);
template = template.*(K./sum(sum(template)));
CIJ = zeros(N);

while ((sum(sum(CIJ)) ~= K))
   CIJ = (rand(N)<template);
end;
