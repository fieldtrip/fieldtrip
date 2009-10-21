function [i] = nearest(array, val)

% NEAREST return the index of an array nearest to a scalar
% 
% [indx] = nearest(array, val)

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: nearest.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.3  2004/12/06 12:55:57  roboos
% added support for -inf and inf, respectively returning the first and last occurence of the nearest element
%
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights
%

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

