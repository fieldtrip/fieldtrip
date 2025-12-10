function d = pdist2(x, y)

% PDIST2 computes the distance between all points in input matrix X and Y.
%
% Use as
%   d = pdist(x, y)
%
% The output is a matrix with size MX x MY.

% This function was written to be a plugin replacement of the
% similarly-named function in the Matlab stats toolbox.

mx = size(x,1);
my = size(y,1);
d = zeros(mx, my);
for i=1:mx
  for j=1:my
    % compute the Euclidian distance
    d(i,j) = sqrt(sum((x(i,:)-y(j,:)).^2));
  end
end
