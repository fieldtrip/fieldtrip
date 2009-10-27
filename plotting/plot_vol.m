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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

keyvalcheck(varargin, 'forbidden', {'faces', 'edges', 'vertices'});
% get the optional input arguments

faceindex   = keyval('faceindex',   varargin);   if isempty(faceindex),faceindex = 'none';end
vertexindex = keyval('vertexindex',   varargin); if isempty(faceindex),vertexindex ='none';end
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
switch voltype(vol)
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

