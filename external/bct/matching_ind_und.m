function M0 = matching_ind_und(CIJ)
%MATCHING_IND_UND       matching index
%
%   M0 = MATCHING_IND_UND(CIJ) computes matching index for undirected
%   graph specified by adjacency matrix CIJ. Matching index is a measure of
%   similarity between two nodes' connectivity profiles (excluding their
%   mutual connection, should it exist).
%
%   Inputs:     CIJ,    undirected adjacency matrix
%
%   Outputs:    M0,     matching index matrix.
%
%   Richard Betzel, Indiana University, 2013
%
CIJ0 = CIJ;
K = sum(CIJ0);
R = K ~= 0;
N = sum(R);
CIJ = CIJ0(R,R);
I = ~eye(N);
M = zeros(N,N);
for i = 1:N
    
    c1 = CIJ(i,:);
    use = bsxfun(@or,c1,CIJ);
    use(:,i) = 0;
    use = use.*I;
    
    ncon1 = bsxfun(@times,use,c1);
    ncon2 = bsxfun(@times,use,CIJ);
    ncon = sum(ncon1 + ncon2,2);
    
    M(:,i) = 2*sum(ncon1 & ncon2,2)./ncon;
    
end
M = M.*I;
M(isnan(M)) = 0;
M0 = zeros(size(CIJ0));
M0(R,R) = M;