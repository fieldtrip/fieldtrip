function [gx, hx] = splint_gh(x)

% SPLINT_GH implements equations (3) and (5b) of Perrin 1989
% for simultaneous computation of multiple values
%
% Use as
%   [gx, hx] = splint_gh(x)

% Copyright (C) 2004-2009, Robert Oostenveld
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

M = 4;          % constant in denominator
N = 9;          % number of terms for series expansion
p  = zeros(1,N);
gx = zeros(size(x));
hx = zeros(size(x));
x(find(x>1)) = 1;       % to avoid rounding off errors
x(find(x<-1)) = -1;     % to avoid rounding off errors
for i=1:size(x,1)
  for j=1:size(x,2)
    for k=1:N
      p(k) = plgndr(k,0,x(i,j));
    end
    gx(i,j) =  sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^M     .* p) / (4*pi);
    hx(i,j) = -sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^(M-1) .* p) / (4*pi);
  end
end

