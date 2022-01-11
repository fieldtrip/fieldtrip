% w = hann(n)
%   see hanning: it's almost the same, apart from the fact that a hann has
%   a 0 at the edge, and a symmetric hanning doesn't
function w = hann(n, opt)

if nargin<2
  opt = 'symmetric';
end

if isequal(opt, 'symmetric')
  w = hanning(n-2, opt);
  w = [0;w;0];
elseif isequal(opt, 'periodic')
  w = hanning(n, opt);
end

% overrule the previous output
if n==1
  w = 1;
end
