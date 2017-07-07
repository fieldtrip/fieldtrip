function s = pad(s, n, side, c)

% PAD adds leading or trailing characters, such as spaces, to the left or
% right of an existing string.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

if nargin<4 || isempty(c)
  c = ' ';
end

if nargin<3 || isempty(side)
  side = 'right';
end

if nargin<2 || isempty(n)
  if iscell(s)
    n = max(cellfun(@length, s(:)));
  else
    n = length(s);
  end
end

if iscell(s)
  % use recursion to deal with cell-arrays
  for i=1:numel(s)
    s{i} = pad(s{i}, n, side, c);
  end
  
else
  % this is where the actual work happens
  assert(size(s,1)<2);
  assert(ischar(s));
  assert(numel(c)==1);
  assert(ischar(c));
  
  if length(s)>=n
    return
  else
    c = repmat(c, 1, n-length(s));
    switch (side)
      case 'left'
        s = [c s];
      case 'right'
        s = [s c];
      otherwise
        error('unsupported side')
    end % switch
  end
end
