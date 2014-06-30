function bool = inside_contour(pos, contour)

npos = size(pos,1);
ncnt = size(contour,1);
x = pos(:,1);
y = pos(:,2);

minx = min(x);
miny = min(y);
maxx = max(x);
maxy = max(y);

bool = true(npos,1);
bool(x<minx) = false;
bool(y<miny) = false;
bool(x>maxx) = false;
bool(y>maxy) = false;

% the summed angle over the contour is zero if the point is outside, and 2*pi if the point is inside the contour
% leave some room for inaccurate f
critval = 0.1;

% the remaining points have to be investigated with more attention
sel = find(bool);
for i=1:length(sel)
  contourx = contour(:,1) - pos(sel(i),1);
  contoury = contour(:,2) - pos(sel(i),2);
  angle = atan2(contoury, contourx);
  % angle = unwrap(angle);
  angle = my_unwrap(angle);
  total = sum(diff(angle));
  bool(sel(i)) = (abs(total)>critval);
end

function x = my_unwrap(x)
% this is a faster implementation of the MATLAB unwrap function
% with hopefully the same functionality
d    = diff(x);
indx = find(abs(d)>pi);
for i=indx(:)'
  if d(i)>0
    x((i+1):end) = x((i+1):end) - 2*pi;
  else
    x((i+1):end) = x((i+1):end) + 2*pi;
  end
end


