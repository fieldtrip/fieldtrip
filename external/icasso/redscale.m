function map=redscale(N,B)
%function map=redscale(N,B)
%
%PURPOSE
%
%Makes a colormap of shades of red from white to red.
%
%INPUT
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% N   (scalar) number of colors
% [B] (scalar) the Red value (0...1)
%      to start from default: 1 (bright red) 
%
%OUTPUT
%
%map (matrix) a Nx3 matrix of RGB values
%

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 johan 100105

% scale factor
F=0.5;

if nargin<2|isempty(B),
  B=1;
end

% define colors

map(:,1)=linspace(1,B.^2,N)';
map(:,2)=linspace(1,0,N)';
map(:,3)=map(:,2);
map(:,1)=map(:,1).^F;
map(:,2:3)=map(:,2:3).^F;
