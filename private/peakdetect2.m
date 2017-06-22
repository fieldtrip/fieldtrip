function [p, v] = peakdetect2(dat, val, mindist)

% PEAKDETECT2 detects peaks above a certain threshold in single-channel data
%
% Use as
%   [pindx, pval] = peakdetect(signal, min, mindist)
%
% mindist is optional, default is 1
%
% See also PEAKDETECT, PEAKDETECT3

% Copyright (C) 2000, Robert Oostenveld
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

if nargin<3
  mindist=1;
end

i = find(dat>val);
m = dat(i);
d = diff(i);
jump = (d>mindist);
p = [];

sect=1;
while sect<=length(d)
  if jump(sect)
    p = [p i(sect)];
  else
    s = [];
    while ~jump(sect) && sect<length(d)
      s = [s sect];
      sect = sect + 1;
    end
    [lm, li] = max(m(s));
    p = [p i(s(li))];
  end
  sect = sect+1;
end

if nargout>1
  v = dat(p);
end

