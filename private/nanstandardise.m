function [x,mx,sx] = standardise(x,dim)

% X = NANSTANDARDISE(X, DIM) computes the zscore of a matrix along dimension 
% dim, taking nans into account
% has similar functionality as the stats-toolbox's zscore function

% Copyright (C) 2010, Jan-Mathijs Schoffelen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin == 1, dim = find(size(x)>1,1,'first'); end

siz    = size(x);
n      = sum(~isnan(x),dim);
x(isnan(x)) = 0;
repsiz = siz;
repsiz(setdiff(1:numel(siz), dim)) = 1;
ressiz      = [siz 1];
if dim>1,
  ressiz(dim) = [];
else
  ressiz(dim) = 1;
end
mx     = sum(x,dim)./n;
x      = x - repmat(mx, repsiz);
%mx     = reshape(mx, ressiz);
sx     = sqrt(sum(x.^2,dim)./n);
x      = x ./repmat(sx, repsiz);
%sx     = reshape(sx, ressiz);
