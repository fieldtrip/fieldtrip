function [d] = dist(x);

% DIST computes the euclidian distance between the columns of the input matrix
%
% Use as
%   [d] = dist(x)
% where x is for example Nx3 for points in 3D space.
%
% This function serves as a replacement for the dist function in the Neural
% Networks toolbox.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: dist.m,v $
% Revision 1.1  2005/05/26 07:31:54  roboos
% new implementation to replace Neural Networks version (required by createlayout)
%

n = size(x,2);
d = zeros(n,n);
for i=1:n
  for j=(i+1):n
    d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
    d(j,i) = d(i,j);
  end
end

