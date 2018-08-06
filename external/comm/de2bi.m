function b = de2bi(d, n, p, flg)

% DE2BI converts decimal numbers to their binary representation
%
% This is a drop-in replacement for the corresponding function from
% the Mathworks Communications System Toolbox.

% see http://se.mathworks.com/help/releases/R2015b/comm/ref/de2bi.html

d = d(:);
m = numel(d);

if nargin<2 || isempty(n)
  n = floor(log(max(d))/log(2)+1);
end

if nargin>=3 && p~=2
  % only supported for base 2
  error('not implemented');
end

if nargin<4 || isempty(flg)
  flg = 'right-msb';
end

b = zeros(m,n);
for i=1:m
  b(i,:) = bitget(d(i),1:n);
end

switch flg
  case 'left-msb'
    % flip left-right
    b = fliplr(b);
  case 'right-msb'
    % keep as it is
  otherwise
    error('not implemented');
end

