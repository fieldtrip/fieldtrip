function [InvRtrunc,q,u,v]=nmt_svdtruncinv(R,signalspace)
% [Rtrunc,q,u,v] = nmt_svdtruncinv(R,signalspace)
%
% Allows user to define signalspace of an arbitrary matrix R and perform an
% SVD-based pseudoinverse using knowledge of the defined signalspace.
% This provides more control than the pinv function and often better
% results than Tikhonov regularization for heavily rank-deficient matrices.
%
% signalspace should be vector of which SVD components to include, e.g. [1:90]
% to set rank to 90.

[u,q,v]=svd(R,'econ');
qinv = diag(1./diag(q));
notsignalspace = setdiff(1:length(q),signalspace);
qinv(notsignalspace,notsignalspace)=0;
InvRtrunc = v*qinv*u';
