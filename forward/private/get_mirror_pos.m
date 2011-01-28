function P2 = get_mirror_pos(P1,vol);

% GET_MIRROR_POS calculates the position of a point symmetric to pnt with respect to a plane
%       
% P2 = get_mirror_pos(P1,vol);

% Copyright (C) 2011, Cristiano Micheli 
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
% $Id: get_mirror_pos.m $

P2 = [];

% define the plane
pnt = vol.pnt;
ori = vol.ori;
% plane = defplane(pnt,ori);
% P2 = symm(plane,P1);
