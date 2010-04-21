function plot_vol(vol, varargin)

% PLOT_VOL visualizes the boundaries in the vol structure constituting the
% geometrical information of the forward model
%
% Use as
%   hs = plot_vol(vol, varargin)
%
% Graphic facilities are available for vertices, edges and faces. A list of
% the arguments is given below with the correspondent admitted choices.
%
%     'facecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'vertexcolor'   [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'edgecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'faceindex'     true or false
%     'vertexindex'   true or false
%
% Example
%   vol.r = [86 88 92 100];
%   vol.o = [0 0 40];
%   figure, plot_vol(vol)

% Copyright (C) 2009, Cristiano Micheli
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

warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
keyvalcheck(varargin, 'forbidden', {'faces', 'edges', 'vertices'});
faceindex   = keyval('faceindex',   varargin);   if isempty(faceindex),faceindex = 'none';end
vertexindex = keyval('vertexindex',   varargin); if isempty(vertexindex),vertexindex ='none';end
vertexsize  = keyval('vertexsize',    varargin); if isempty(vertexsize),  vertexsize = 10;    end
facecolor   = keyval('facecolor',     varargin); if isempty(facecolor),facecolor = 'white'; end 
vertexcolor = keyval('vertexcolor',   varargin); if isempty(vertexcolor),vertexcolor ='none';end
edgecolor   = keyval('edgecolor',     varargin); if isempty(edgecolor),edgecolor = 'k';end
facealpha   = keyval('facealpha',     varargin); if isempty(facealpha),facealpha = 1;end 
map         = keyval('colormap',      varargin);

faceindex   = istrue(faceindex);
vertexindex = istrue(vertexindex);

% we will probably need a sphere, so let's prepare one
[pnt, tri] = icosahedron162;

% prepare a single or multiple triangulated boundaries
switch ft_voltype(vol)
  case {'singlesphere' 'concentric'}
    vol.r = sort(vol.r);
    bnd = [];
    for i=1:length(vol.r)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(3);
      bnd(i).tri = tri;
    end

  case 'multisphere'
    bnd = [];
    for i=1:length(vol.label)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(i,1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(i,2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(i,3);
      bnd(i).tri = tri;
    end

  case {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'nolte'}
    % these already contain one or multiple triangulated surfaces for the boundaries
    bnd = vol.bnd;

  otherwise
    error('unsupported voltype')
end

 
% plot the triangulated surfaces of the volume conduction model
for i=1:length(bnd)
  plot_mesh(bnd(i),'faceindex',faceindex,'vertexindex',vertexindex, ...
    'vertexsize',vertexsize,'facecolor',facecolor,'edgecolor',edgecolor, ...
    'vertexcolor',vertexcolor,'facealpha',facealpha);
end

