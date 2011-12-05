function C=clustering_coef_wd(W)
%CLUSTERING_COEF_WD     Clustering coefficient
%
%   C = clustering_coef_wd(W);
%
%   The weighted clustering coefficient is the average "intensity" of 
%   triangles around a node.
%
%   Input:      W,      weighted directed connection matrix
%
%   Output:     C,      clustering coefficient vector
%
%   Reference: Fagiolo (2007) Phys Rev E 76:026107.
%
%
%   Mika Rubinov, UNSW, 2007-2010


%Methodological note (also see clustering_coef_bd)
%The weighted modification is as follows:
%- The numerator: adjacency matrix is replaced with weights matrix ^ 1/3
%- The denominator: no changes from the binary version
%
%The above reduces to symmetric and/or binary versions of the clustering 
%coefficient for respective graphs.

A=W~=0;                     %adjacency matrix
S=W.^(1/3)+(W.').^(1/3);	%symmetrized weights matrix ^1/3
K=sum(A+A.',2);            	%total degree (in + out)
cyc3=diag(S^3)/2;           %number of 3-cycles (ie. directed triangles)
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
CYC3=K.*(K-1)-2*diag(A^2);	%number of all possible 3-cycles
C=cyc3./CYC3;               %clustering coefficient

