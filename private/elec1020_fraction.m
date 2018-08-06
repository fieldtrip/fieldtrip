function [pnt] = elec1020_fraction(cnt1, cnt2, fraction)

% ELEC1020_FRACTION

% Copyright (C) 2003, Robert Oostenveld
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

ncnt = size(cnt1,1);

% determine the total length of the contour
tot_l = 0;
for i=1:ncnt
  tot_l = tot_l + pntdist(cnt1(i,:), cnt2(i,:));
end

frac_l = fraction * tot_l;

% propagate along the contour untill we get at the desired fraction
sum_l = 0;
for i=1:ncnt
  seg_l = pntdist(cnt1(i,:), cnt2(i,:));
  if (sum_l+seg_l)>=frac_l
    % the desired point lies on this segment
    la = frac_l - sum_l;
    vec = cnt2(i,:)-cnt1(i,:);
    pnt = cnt1(i,:) + la * vec/norm(vec);
    return
  else
    sum_l = sum_l + seg_l;
    sum_f = sum_l/tot_l;
  end
end

pnt = [nan nan nan];

