function test_nearest

% TEST test_nearest
% TEST nearest

% Use as
%   [indx] = nearest(array, val, insideflag, toleranceflag)

% these are some normal cases
assert(nearest([1 2 3], 1)==1)
assert(nearest([1 2 3], 2)==2)
assert(nearest([1 2 3], 3)==3)
assert(nearest([1 2 3], -inf)==1)
assert(nearest([1 2 3],  inf)==3)

% unsorted arrays should be supported
assert(nearest([1 3 2], -inf)==1)
assert(nearest([3 1 2], -inf)==2)
assert(nearest([3 1 2], 1)==2)

% outside the range
try
  nearest([1 2 3], 0, true);
  error('this should have returned an error');
end
try
  nearest([1 2 3], 4, true);
  error('this should have returned an error');
end

% just inside the tolerance range
assert(nearest([1 2 3], 0.5, true, true)==1);
assert(nearest([1 2 3], 3.5, true, true)==3);

% just outside the tolerance range
try
  nearest([1 2 3], 0.499, true, true);
  error('this should have returned an error');
end
try
  nearest([1 2 3], 3.501, true, true);
  error('this should have returned an error');
end

