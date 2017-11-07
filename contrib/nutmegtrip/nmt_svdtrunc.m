function [Rtrunc,q,u,v]=nmt_svdtrunc(R,signalspace)
% [Rtrunc,q,u,v] = nmt_svdtrunc(R,signalspace)
%
% Reduces the rank of an arbitrary matrix. 
% signalspace should be vector of which components to include, e.g. [1 2]
% for the equivalent of cfg.reducerank=2

[u,q,v]=svd(R,'econ');

notsignalspace = setdiff(1:length(q),signalspace);
q(notsignalspace,notsignalspace)=0;
Rtrunc = u*q*v';
