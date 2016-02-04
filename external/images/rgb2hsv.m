function h = rgb2hsv(m)

% RGB2HSV converts red-green-blue colors to hue-saturation-value.
%
% this code is based on the comments in
% http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both

if isvector(m) && length(m)==3
  r = m(1);
  g = m(2);
  b = m(3);
elseif ismatrix(m) && size(m,2)==3
  r = m(:,1);
  g = m(:,2);
  b = m(:,3);
elseif ndims(m)==3
  r = m(:,:,1);
  g = m(:,:,2);
  b = m(:,:,3);
end

% ensure these to be double precision column vectors
r = double(r(:));
g = double(g(:));
b = double(b(:));

minrgb = min([r g b], [], 2);
maxrgb = max([r g b], [], 2);
delta  = maxrgb - minrgb;

h = nan(size(delta)); % will be determined further down
v = maxrgb;
s = delta./maxrgb;

% set for each element in the array the proper hue, depending on whether red, green or blue is the largest value
% the 0, 2 and 4 are there to end up in the right corner of the color circle
sel = (r==maxrgb);
hue = 0.0 + (g-b)./delta; % between yellow & magenta
h(sel) = hue(sel);
sel = (g==maxrgb);
hue = 2.0 + (b-r)./delta; % between cyan & yellow
h(sel) = hue(sel);
sel = (b==maxrgb);
hue = 4.0 + (r-g)./delta; % between magenta & cyan
h(sel) = hue(sel);

% these should not be NaN
h(delta==0) = 0;
s(delta==0) = 0;

% convert to degrees between 0 and 360
h = h * 60;
h = mod(h, 360);

% convert to a fraction between 0 and 1
h = h/360;

if nargout==3
  % return the three output arguments separately
else
  % combine the hue-saturation-value in a single array
  h = reshape([h(:) s(:) v(:)], size(m));
  clear s v
end
