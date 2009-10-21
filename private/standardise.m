function [x,mx,sx] = standardise(x,dim)

% X = STANDARDISE(X, DIM) computes the zscore of a matrix along dimension dim
% has extended functionality as compared to the stats-toolbox's zscore function

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: standardise.m,v $
% Revision 1.2  2009/06/16 15:44:29  jansch
% added automatic dim detection (first non singular dimension) for nargin==1.
% added mean and std to output
%
% Revision 1.1  2009/05/19 15:59:11  jansch
% first commitment into cvs
%

if nargin == 1, dim = find(size(x)>1,1,'first'); end

n      = size(x,dim);
mx     = mean(x,dim);
%sx     = std(x,0,dim);
sx     = std(x,1,dim);
repvec = ones(1,length(size(x)));
repvec(dim) = n; 
x      = (x - repmat(mx,repvec))./repmat(sx,repvec);
