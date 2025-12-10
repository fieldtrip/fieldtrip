function test_nearest

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY nearest
% DATA no

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


try
  nearest(.1:.1:1.0,[3 8]);
catch me
  if ~strcmp(me.message,'The selected value should be a scalar, non-scalar value is not supported anymore due to ambiguity in functional behavior')
    error('wrong error message in nearest')
  end
end

% create a large array and test

x = 0:492:(492*(25*10^6));
y = 8424306*492;
indx = nearest(x,y);
indx2 = find(x<=y, 1, 'last');

% temporal comparison

tic
for i=1:10
  indx = nearest(x,y);
end
t1 = toc;

tic
for i=1:10
  x = sort(x);
  indx2 = find(x<=y, 1, 'last');
  indx3 = find(x>=y, 1, 'first');
  if abs(x(indx2)-y) < abs(x(indx3)-y)
    indx4 = indx2;
  else
    indx4 = indx3;
  end
end
t2 = toc;

fprintf('Time needed for nearest function: %.2g s per 10 calls\n', t1/10);
fprintf('Time needed for find function: %.2g s per 10 calls\n', t2/10);

if indx~=indx2
  error('nearest does not output the correct value');
  fprintf('nearest is off by %d samples', indx2-indx);
end
% indx2 = 8424307



