function [indx] = nearest(array, val, insideflag, toleranceflag)

% NEAREST return the index of an array nearest to a scalar
%
% Use as
%   [indx] = nearest(array, val, insideflag, toleranceflag)
%
% The second input val can be a scalar, or a [minval maxval] vector for
% limits selection.
%
% If not specified or if left empty, the insideflag and the toleranceflag
% will default to false.
%
% The boolean insideflag can be used to specify whether the value should be
% within the array or not. For example nearest(1:10, -inf) will return 1,
% but nearest(1:10, -inf, true) will return an error because -inf is not
% within the array.
%
% The boolean toleranceflag is used when insideflag is true. It can be used
% to specify whether some tolerance should be allowed for values that are
% just outside the array. For example nearest(1:10, 0.99, true, false) will
% return an error, but nearest(1:10, 0.99, true, true) will return 1. The
% tolerance that is allowed is half the distance between the subsequent
% values in the array.

% Copyright (C) 2002-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

mbreal(array);
mbreal(val);

mbvector(array);
assert(all(~isnan(val)), 'incorrect value (NaN)');

if numel(val)==2
  % interpret this as a range specification like [minval maxval]
  % see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1431
  intervaltol=eps;
  sel = find(array>=val(1) & array<=val(2));
  if isempty(sel)
    error('The limits you selected are outside the range available in the data');
  end
  indx(1) = sel(1);
  indx(2) = sel(end);
  if indx(1)>1 && abs(array(indx(1)-1)-val(1))<=intervaltol
    indx(1)=indx(1)-1;
  end
  if indx(2)<length(array) && abs(array(indx(2)+1)-val(2))<=intervaltol
    indx(2)=indx(2)+1;
  end
  return
end

mbscalar(val);

if nargin<3 || isempty(insideflag)
  insideflag = false;
end

if nargin<4 || isempty(toleranceflag)
  toleranceflag = false;
end

% ensure that it is a column vector
array = array(:);

% determine the most extreme values in the array
minarray = min(array);
maxarray = max(array);

% do some strict checks whether the value lies within the min-max range
if insideflag
  if ~toleranceflag
    if val<minarray || val>maxarray
      if numel(array)==1
        warning('the selected value %g should be within the range of the array from %g to %g', val, minarray, maxarray);
      else
        error('the selected value %g should be within the range of the array from %g to %g', val, minarray, maxarray);
      end
    end
  else
    if ~isequal(array, sort(array))
      error('the input array should be sorted from small to large');
    end
    if numel(array)<2
      error('the input array must have multiple elements to compute the tolerance');
    end
    mintolerance = (array(2)-array(1))/2;
    maxtolerance = (array(end)-array(end-1))/2;
    if val<(minarray-mintolerance) || val>(maxarray+maxtolerance)
      error('the value %g should be within the range of the array from %g to %g with a tolerance of %g and %g on both sides', val, minarray, maxarray, mintolerance, maxtolerance);
    end
  end % toleragceflag
end % insideflag

% FIXME it would be possible to do some soft checks and potentially give a
% warning in case the user did not explicitly specify the inside and
% tolerance flags

% note that [dum, indx] = min([1 1 2]) will return indx=1
% and that  [dum, indx] = max([1 2 2]) will return indx=2
% whereas it is desired to have consistently the match that is most towards the side of the array

if val>maxarray
  % return the last occurence of the largest number
  [dum, indx] = max(flipud(array));
  indx = numel(array) + 1 - indx;
  
elseif val<minarray
  % return the first occurence of the smallest number
  [dum, indx] = min(array);
  
else
  % implements a threshold to correct for errors due to numerical precision
  % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=498 and http://bugzilla.fcdonders.nl/show_bug.cgi?id=1943
  %   if maxarray==minarray
  %     precision = 1;
  %   else
  %     precision = (maxarray-minarray) / 10^6;
  %   end
  %
  %   % return the first occurence of the nearest number
  %   [dum, indx] = min(round((abs(array(:) - val)./precision)).*precision);
  
  % use find instead, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1943
  wassorted = true;
  if ~issorted(array)
    wassorted = false;
    [array, xidx] = sort(array);
  end
  
  indx2 = find(array<=val, 1, 'last');
  indx3 = find(array>=val, 1, 'first');
  if abs(array(indx2)-val) <= abs(array(indx3)-val)
    indx = indx2;
  else
    indx = indx3;
  end
  if ~wassorted
    indx = xidx(indx);
  end
  
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
if ndims(a) > 2 || (size(a, 1) > 1 && size(a, 2) > 1)
  error('Argument to mbvector must be a vector');
end

