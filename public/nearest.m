function [i] = nearest(array, val)

% NEAREST return the index of an array nearest to a scalar
% 
% [indx] = nearest(array, val)

% Copyright (C) 2002, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

mbreal(array);
mbreal(val);

mbvector(array);
mbscalar(val);

% ensure that it is a column vector
array = array(:);

if isnan(val)
  error('incorrect value')
end

if val>max(array)
  % return the last occurence of the nearest number
  [dum, i] = max(flipud(array));
  i = length(array) + 1 - i;
else
  % return the first occurence of the nearest number
  [mindist, i] = min(abs(array(:) - val));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbreal(a)
if ~isreal(a)
  error('Argument to mbreal must be real');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbscalar(a)
if ~all(size(a)==1)
  error('Argument to mbscalar must be scalar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbvector(a)
if ndims(a) > 2 | (size(a, 1) > 1 & size(a, 2) > 1)
  error('Argument to mbvector must be a vector');
end

