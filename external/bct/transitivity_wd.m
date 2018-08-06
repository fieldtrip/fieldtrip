function T=transitivity_wd(W)
%TRANSITIVITY_WD    Transitivity
%
%   T = transitivity_wd(W);
%
%   Transitivity is the ratio of 'triangles to triplets' in the network.
%   (A classical version of the clustering coefficient).
%
%   Input:      W       weighted directed connection matrix
%
%   Output:     T       transitivity scalar
%
%   Note:       All weights must be between 0 and 1.
%               This may be achieved using the weight_conversion.m function,
%               W_nrm = weight_conversion(W, 'normalize');
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
%   2015: Expanded documentation


%   Methodological note (also see note for clustering_coef_bd)
%   The weighted modification is as follows:
%   - The numerator: adjacency matrix is replaced with weights matrix ^ 1/3
%   - The denominator: no changes from the binary version
%
%   The above reduces to symmetric and/or binary versions of the clustering
%   coefficient for respective graphs.

A    = W~=0;                        % adjacency matrix
S    = W.^(1/3)+(W.').^(1/3);      	% symmetrized weights matrix ^1/3
K    = sum(A+A.',2);               	% total degree (in + out)
cyc3 = diag(S^3)/2;                	% number of 3-cycles (ie. directed triangles)
CYC3 = K.*(K-1)-2*diag(A^2);        % number of all possible 3-cycles
T    = sum(cyc3)./sum(CYC3);       	% transitivity

