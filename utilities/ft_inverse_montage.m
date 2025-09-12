function output = ft_inverse_montage(input)

% FT_INVERSE_MONTAGE constructs an inverse of the input montage matrix,
% allowing a linear projection of the channel timeseries to be undone.
%
% Use as
%   montage = ft_inverse_montage(montage)
%
% See also FT_APPLY_MONTAGE

% Copyright (C) 2025, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% ensure up-to-date field names
input = fixoldorg(input);

% swap the role of the old and new channels
output.labelnew    = input.labelold;
output.labelold    = input.labelnew;
if isfield(input, 'chantypeold')
  output.chantypenew = input.chantypeold;
end
if isfield(input, 'chantypenew')
  output.chantypeold = input.chantypenew;
end
if isfield(input, 'chanunitold')
  output.chanunitnew = input.chanunitold;
end
if isfield(input, 'chanunitnew')
  output.chanunitold = input.chanunitnew;
end

% ensure it is not a sparse or a single representation
input.tra = double(full(input.tra));

if rank(input.tra) < length(input.tra)
  ft_warning('the linear projection for the montage is not full-rank, the result may have reduced dimensionality');
  output.tra = pinv(input.tra);
else
  output.tra = inv(input.tra);
end
