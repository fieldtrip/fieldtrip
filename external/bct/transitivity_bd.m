function T=transitivity_bd(A)
%TRANSITIVITY_BD    Transitivity
%
%   T = transitivity_bd(A);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      A       binary directed connection matrix
%
%   Output:     T       transitivity scalar
%
%   Reference:  Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
%               based on Fagiolo (2007) Phys Rev E 76:026107.
%
%
%   Contributors:
%   Mika Rubinov, UNSW/University of Cambridge
%   Christoph Schmidt, Friedrich Schiller University Jena
%   Andrew Zalesky, University of Melbourne
%   2007-2015

%   Modification history:
%   2007: original (MR)
%   2013, 2015: removed tests for absence of nodewise 3-cycles (CS,AZ)

%   Methodological note: In directed graphs, 3 nodes generate up to 8 
%   triangles (2*2*2 edges). The number of existing triangles is the main 
%   diagonal of S^3/2. The number of all (in or out) neighbour pairs is 
%   K(K-1)/2. Each neighbour pair may generate two triangles. "False pairs"
%   are i<->j edge pairs (these do not generate triangles). The number of 
%   false pairs is the main diagonal of A^2. Thus the maximum possible 
%   number of triangles = (2 edges)*([ALL PAIRS] - [FALSE PAIRS])
%                       = 2 * (K(K-1)/2 - diag(A^2))
%                       = K(K-1) - 2(diag(A^2))

S    = A+A.';                           % symmetrized input graph
K    = sum(S,2);                        % total degree (in + out)
cyc3 = diag(S^3)/2;                     % number of 3-cycles (ie. directed triangles)
CYC3 = K.*(K-1)-2*diag(A^2);            % number of all possible 3-cycles
T    = sum(cyc3)./sum(CYC3);            % transitivity
