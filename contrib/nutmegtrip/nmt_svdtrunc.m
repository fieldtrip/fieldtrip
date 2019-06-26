function [Rtrunc,q,u,v]=nmt_svdtrunc(R,signalspace)
% [Rtrunc,q,u,v] = nmt_svdtrunc(R,signalspace)
%
% Allows user to reject undesired SVD components of an arbitrary matrix R.
% It can be used to, e.g., reduce the matrix's rank.
%
% signalspace should be vector of which components to include, e.g. [1 2]
% for the equivalent of cfg.reducerank=2

[u,q,v]=svd(R,'econ');

notsignalspace = setdiff(1:length(q),signalspace);
q(notsignalspace,notsignalspace)=0;
Rtrunc = u*q*v';
