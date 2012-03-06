function [dv]=dvComp(Xtst,Xtrn,kernel,alphab,varargin)
% Compute the decision values for a kernel classifier
%
% [dv]=dvComp(Xtst,kernel,Xtrn,alphab,...)
%
% dv = K(Xtrn,Xtst)*alphab(1:end-1)+alphab(end)
%
% Inputs:
%  Xtst   -- [N x d] test set
%  Xtrn   -- [Ntrn x d] training set
%  kernel -- the kernel function *as used in training*, 
%            will be passed to compKernel
%  alphab -- [Ntrn+1 x 1] set of decision function parameters [alpha;b]
%  ...    -- additional parameters to pass to compKernel
%
% Outputs:
%  dv     -- [N x 1] decision values
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied

K  = compKernel(Xtst,Xtrn,kernel,varargin{:});
dv = K*alphab(1:end-1)+alphab(end);
return;

%----------------------------------------------------------
function testCase()
[X,Y]=mkMultiClassTst([-1 0; 1 0; .2 .5],[400 400 50],[.3 .3; .3 .3; .2 .2],[],[-1 1 1]);[dim,N]=size(X);

trnInd = randn(N,1)>0;
K=compKernel(X(:,trnInd)',[],'nlinear');
[alphab,f,J]=rkls(K,Y(trnInd),1);

dv = dvComp(X(:,trnInd)',[],'nlinear',alphab);
max(abs(dv-f))

dv = dvComp(X(:,~trnInd)',X(:,trnInd)','nlinear',alphab);