function [x,mx,sx] = standardise(x,dim)

% X = STANDARDISE(X, DIM) computes the zscore of a matrix along dimension dim
% has extended functionality as compared to the stats-toolbox's zscore function

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin == 1, dim = find(size(x)>1,1,'first'); end

n      = size(x,dim);
mx     = mean(x,dim);
%sx     = std(x,0,dim);
sx     = std(x,1,dim);
repvec = ones(1,length(size(x)));
repvec(dim) = n; 
x      = (x - repmat(mx,repvec))./repmat(sx,repvec);
