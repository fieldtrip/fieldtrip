function [cfg] = ft_headmodelplot(cfg, data)

% FT_HEADMODELPLOT makes a 3D visualisation of the volume conductor model
% and optionally of the gradiometer positions and headshape. It can
% be used for example to check CTF multiple-sphere head models.
%
% Use as
%   ft_headmodelplot(cfg)
%   ft_headmodelplot(cfg, data)
%
% You should specify the volume conductor model with
%   cfg.hdmfile       = string, file containing the volume conduction model
% or alternatively
%   cfg.vol           = structure with volume conduction model
%
% If the sensor information is obtained by FT_FETCH_SENS.
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = Nx3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
% You can also use the FT_PREPARE_LEADFIELD or FT_SOURCEANALYSIS functions
% to create a grid with dipole positions.
%
% Other options are
%   cfg.channel          = cell-array, see FT_CHANNELSELECTION
%   cfg.spheremesh       = number of vertices for spheres, either 42, 162 or 642
%   cfg.plotheadsurface  = 'yes' or 'no', is constructed from head model
%   cfg.plotbnd          = 'yes' or 'no'
%   cfg.plotspheres      = 'yes' or 'no'
%   cfg.plotspherecenter = 'yes' or 'no'
%   cfg.plotgrid         = 'yes' or 'no'
%   cfg.plotinside       = 'yes' or 'no'
%   cfg.plotoutside      = 'yes' or 'no'
%   cfg.plotsensors      = 'yes' or 'no' plot electrodes or gradiometers
%   cfg.plotcoil         = 'yes' or 'no' plot all gradiometer coils
%   cfg.plotlines        = 'yes' or 'no' plot lines from sensor to head surface
%   cfg.surftype         = 'edges'or 'faces'
%   cfg.headshape        = a filename containing headshape, a structure containing a
%                          single triangulated boundary, or a Nx3 matrix with surface
%                          points
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_VOLUMEDOWNSAMPLE FT_VOLUMENORMALISE FT_VOLUMERESLICE
% FT_VOLUMEWRITE FT_VOLUMELOOKUP FT_VOLUMEREALIGN, FT_VOLUMESEGMENT,
% FT_SOURCEANALYSIS

% Undocumented local options:
% cfg.surface_facecolor
% cfg.surface_edgecolor
% cfg.surface_facealpha

%
% This function depends on FT_PREPARE_VOL_SENS which has the following options:
% cfg.channel, documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, (default set in FT_HEADMODELPLOT at cfg.vol =[]), documented

% Copyright (C) 2004-2007, Robert Oostenveld
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

revision = '$Id$';

% throw a warning about this function being deprecated
warning('The ft_headmodelplot function is deprecated and is not supported anymore by the FieldTrip team. Please look at this FAQ http://fieldtrip.fcdonders.nl/faq/how_can_i_visualize_the_different_geometrical_objects_that_are_needed_for_forward_and_inverse_computations for the suggested alternative approach. If you still want to use this function and get errors, you can manually move this function from its present location into the fieldtrip main directory, where it still may work');

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% these are suitable RGB colors
skin_surface   = [255 213 119]/255;
outer_skull_surface  = [140  85  85]/255;
inner_skull_surface  = [202 100 100]/255;
cortex = [255 213 119]/255;

% set the defaults
if ~isfield(cfg, 'surface_facecolor'), cfg.surface_facecolor = skin_surface;   end
if ~isfield(cfg, 'surface_edgecolor'), cfg.surface_edgecolor = 'none';         end
if ~isfield(cfg, 'surface_facealpha'), cfg.surface_facealpha = 0.7;            end
if ~isfield(cfg, 'surftype'),          cfg.surftype = 'faces';                 end

if ~isfield('data', 'var')    
  % this will be passed into the prepare_headmodel function further down
  data = [];
end

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});

if ~isfield(cfg, 'vol') && ~isfield(cfg, 'hdmfile')
  cfg.vol = [];  % FIXME why is this empty setting neccessary?
end

% set the defaults that apply both to EEG and MEG
if ~isfield(cfg, 'spheremesh'),       cfg.spheremesh = 642;          end
if ~isfield(cfg, 'plotsensors'),      cfg.plotsensors = 'yes';       end
if ~isfield(cfg, 'plotheadsurface'),  cfg.plotheadsurface = 'yes';   end
if ~isfield(cfg, 'plotgrid'),         cfg.plotgrid = 'yes';          end
if ~isfield(cfg, 'plotinside'),       cfg.plotinside = 'yes';        end
if ~isfield(cfg, 'plotoutside'),      cfg.plotoutside = 'no';        end
if ~isfield(cfg, 'plotbnd'),          cfg.plotbnd = 'no';            end
if ~isfield(cfg, 'plotfiducial'),     cfg.plotfiducial = 'no';       end

% extract/read the gradiometer and volume conductor
[vol, sens, cfg] = prepare_headmodel(cfg, data);

if strcmp(cfg.plotgrid, 'yes')
  if isfield(cfg.grid, 'pos')
    % use the specified grid
    sourcegrid.pos      = cfg.grid.pos;
    sourcegrid.inside   = cfg.grid.inside;
    sourcegrid.outside  = cfg.grid.outside;
  else
    % construct the dipole grid according to the configuration
    tmpcfg = [];
    tmpcfg.vol  = vol;
    tmpcfg.grad = sens; % this can be electrodes or gradiometers
    % copy all options that are potentially used in ft_prepare_sourcemodel
    try, tmpcfg.grid        = cfg.grid;         end
    try, tmpcfg.mri         = cfg.mri;          end
    try, tmpcfg.headshape   = cfg.headshape;    end
    try, tmpcfg.tightgrid   = cfg.tightgrid;    end
    try, tmpcfg.symmetry    = cfg.symmetry;     end
    try, tmpcfg.smooth      = cfg.smooth;       end
    try, tmpcfg.threshold   = cfg.threshold;    end
    try, tmpcfg.spheremesh  = cfg.spheremesh;   end
    try, tmpcfg.inwardshift = cfg.inwardshift;  end
    try, tmpcfg.sourceunits = cfg.sourceunits;  end
    [sourcegrid, tmpcfg] = ft_prepare_sourcemodel(tmpcfg);
  end
else
  % construct an empty dipole grid
  sourcegrid     = [];
  sourcegrid.pos = zeros(0,3);
  sourcegrid.inside = [];
  sourcegrid.outside = [];
end

% determine the type of input data
ismeg          = ft_senstype(sens, 'meg');
iseeg          = ft_senstype(sens, 'eeg');
isbem          = isfield(vol, 'bnd');
issphere       = isfield(vol, 'r');
ismultisphere  = isfield(vol, 'r') && length(vol.r)>4;
issinglesphere = isfield(vol, 'r') && length(vol.r)==1;
isconcentric   = isfield(vol, 'r') && length(vol.r)<=4 && ~issinglesphere;

if ismeg
  % sensors describe MEG data, set the corresponding defaults
  if ~isfield(cfg, 'plotspheres'),      cfg.plotspheres = 'no';        end
  if ~isfield(cfg, 'plotspherecenter'), cfg.plotspherecenter = 'no';   end
  if ~isfield(cfg, 'plotlines'),        cfg.plotlines = 'yes';         end
  if ~isfield(cfg, 'plotcoil'),         cfg.plotcoil = 'no';           end
  if ~isfield(cfg, 'headshape'),        cfg.headshape = [];            end
elseif iseeg
  % sensors describe EEG data, set the corresponding defaults
  if ~isfield(cfg, 'plotspheres'),      cfg.plotspheres = 'no';        end
  if ~isfield(cfg, 'plotspherecenter'), cfg.plotspherecenter = 'no';   end
  if ~isfield(cfg, 'plotlines'),        cfg.plotlines = 'yes';         end
end

chan       = [];
chan.pnt   = sens.chanpos;
chan.label = sens.label;

if issphere
  % determine the number of spheres in the volume model
  Nspheres = length(vol.r);
  fprintf('spherical head model with %d spheres\n', Nspheres);
end

clf
hold on
axis equal
axis vis3d
if ismeg
  axis off
else
  axis on
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
end

if iseeg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotting for EEG
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(cfg.plotgrid, 'yes')
    if strcmp(cfg.plotinside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.inside,1), sourcegrid.pos(sourcegrid.inside,2), sourcegrid.pos(sourcegrid.inside,3), 'k.');
    end
    if strcmp(cfg.plotoutside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.outside,1), sourcegrid.pos(sourcegrid.outside,2), sourcegrid.pos(sourcegrid.outside,3), 'k.');
    end
  end % plotgrid
  
  if strcmp(cfg.plotsensors, 'yes')
    ft_plot_sens(sens, 'style', 'g*');
  end % plotsensors
  
  if strcmp(cfg.plotheadsurface, 'yes')  && ~isempty(vol)
    [pnt, tri] = headsurface(vol, sens);
    h = triplot(pnt, tri, [], cfg.surftype);
    set(h, 'edgecolor', cfg.surface_edgecolor);
    if strcmp(cfg.surftype, 'faces')
      set(h, 'facecolor', cfg.surface_facecolor);
      set(h, 'facealpha', cfg.surface_facealpha);
    end
  end % plotheadsurface
  
  if  strcmp(cfg.plotspheres, 'yes') && ~isempty(vol) && issphere
    % create a triangulated unit sphere
    if cfg.spheremesh==42
      [pnt0, tri] = icosahedron42;
    elseif cfg.spheremesh==162
      [pnt0, tri] = icosahedron162;
    elseif cfg.spheremesh==642
      [pnt0, tri] = icosahedron642;
    end
    Nvertices = size(pnt0,1);
    
    colors = {cortex, inner_skull_surface, outer_skull_surface, skin_surface};
    
    for i=1:Nspheres
      % scale and shift the unit sphere to the proper location
      pnt = pnt0*vol.r(i) + repmat(vol.o, Nvertices,1);
      
      h = triplot(pnt, tri, [], cfg.surftype);
      % set(h, 'FaceVertexCData', 0.5*ones(length(distance),30));
      % set(h, 'FaceVertexCData', [0 0 1]);
      set(h, 'edgecolor', colors{i})
      % set(h, 'edgealpha', 1)
      if strcmp(cfg.surftype, 'faces')
        % set(h, 'AlphaDataMapping', 'direct');
        set(h, 'facealpha', 'interp')
        % set(h, 'facealpha', 0)
        set(h, 'facecolor', 'none');
        % set(h, 'linestyle', 'none');
      end
    end
  end % plotspheres
  
  if  strcmp(cfg.plotbnd, 'yes') && ~isempty(vol) && isbem
    
    Nbnd = numel(vol.bnd);
    
    colors = {skin_surface, outer_skull_surface, inner_skull_surface};
    
    for i=1:Nbnd
      h = triplot(vol.bnd(i).pnt, vol.bnd(i).tri, [], cfg.surftype);
      % set(h, 'FaceVertexCData', 0.5*ones(length(distance),30));
      % set(h, 'FaceVertexCData', [0 0 1]);
      set(h, 'edgecolor', colors{i})
      % set(h, 'edgealpha', 1)
      if strcmp(cfg.surftype, 'faces')
        % set(h, 'AlphaDataMapping', 'direct');
        set(h, 'facealpha', 'interp')
        % set(h, 'facealpha', 0)
        set(h, 'facecolor', 'none');
        % set(h, 'linestyle', 'none');
      end
    end
  end % plotbnd
  
  if strcmp(cfg.plotlines, 'yes') && ~isempty(vol)
    if isbem
      % project the electrodes on the skin surface, on the nearest triangle
      [el] = project_elec(sens.chanpos, vol.bnd(vol.skin_surface).pnt, vol.bnd(vol.skin_surface).tri);
      % this returns [tri, la, mu]
      tri = el(:,1);
      la  = el(:,2);
      mu  = el(:,3);
      for i=1:Nsensors
        v1  = vol.bnd(vol.skin_surface).pnt(vol.bnd(vol.skin_surface).tri(tri(i),1),:);
        v2  = vol.bnd(vol.skin_surface).pnt(vol.bnd(vol.skin_surface).tri(tri(i),2),:);
        v3  = vol.bnd(vol.skin_surface).pnt(vol.bnd(vol.skin_surface).tri(tri(i),3),:);
        prj(i,:) = routlm(v1, v2, v3, la(i), mu(i));
      end
    elseif issphere
      % project the electrodes onto the sphere surface, towards the origin
      nrm = sqrt(sum(sens.chanpos.^2,2));
      prj = vol.r(vol.skin_surface) * sens.chanpos ./ [nrm nrm nrm];
    else
      % in case that no known volume conductor is specified
      prj = sens.chanpos;
    end
    x = [sens.chanpos(:,1) prj(:,1)]';
    y = [sens.chanpos(:,2) prj(:,2)]';
    z = [sens.chanpos(:,3) prj(:,3)]';
    line(x, y, z, 'color', 'm');
  end % plotlines
  
elseif ismeg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plotting for MEG
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(cfg.plotgrid, 'yes')
    if strcmp(cfg.plotinside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.inside,1), sourcegrid.pos(sourcegrid.inside,2), sourcegrid.pos(sourcegrid.inside,3), 'k.');
    end
    if strcmp(cfg.plotoutside, 'yes')
      plot3(sourcegrid.pos(sourcegrid.outside,1), sourcegrid.pos(sourcegrid.outside,2), sourcegrid.pos(sourcegrid.outside,3), 'k.');
    end
  end % plotgrid
  
  if strcmp(cfg.plotsensors, 'yes')
    ft_plot_sens(sens, 'style', 'g*');
  end % plotsensors
  
  if strcmp(cfg.plotcoil, 'yes')
    pnt = sens.coilpos;
    ori = sens.coilori;
    x = [pnt(:,1) pnt(:,1)+ori(:,1)]';
    y = [pnt(:,2) pnt(:,2)+ori(:,2)]';
    z = [pnt(:,3) pnt(:,3)+ori(:,3)]';
    line(x, y, z, 'color', 'g');
    plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'g.');
  end % plotsensors
  
  if ~isempty(cfg.headshape)
    % get the surface describing the head shape
    if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
      % use the headshape surface specified in the configuration
      headshape.pnt = cfg.headshape.pnt;
    elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
      % use the headshape points specified in the configuration
      headshape.pnt = cfg.headshape;
    elseif ischar(cfg.headshape)
      % read the headshape from file
      headshape = ft_read_headshape(cfg.headshape);
    else
      error('cfg.headshape is not specified correctly')
    end
    headshape.pnt = unique(headshape.pnt, 'rows');
    % the triangulation is not used inside this function
    if isfield(headshape, 'tri')
      headshape = rmfield(headshape, 'tri');
    end
    plot3(headshape.pnt(:,1), headshape.pnt(:,2), headshape.pnt(:,3), 'r.');
  end % plotheadshape
  
  if strcmp(cfg.plotspheres, 'yes') && ~isempty(vol) && issphere
    % create a triangulated unit sphere
    if cfg.spheremesh==42
      [pnt0, tri] = icosahedron42;
    elseif cfg.spheremesh==162
      [pnt0, tri] = icosahedron162;
    elseif cfg.spheremesh==642
      [pnt0, tri] = icosahedron642;
    end
    Nvertices = size(pnt0,1);
    for i=1:Nspheres
      % scale and shift the unit sphere to the proper location
      pnt = pnt0*vol.r(i) + repmat(vol.o(i,:),Nvertices,1);
      if Nspheres>4
        distance = sqrt(sum((pnt - repmat(sens.chanpos(i,:),Nvertices,1)).^2, 2));
        distance = (distance-min(distance))/range(distance);
      else
        distance = zeros(Nvertices, 1);
      end
      h = triplot(pnt, tri, [], cfg.surftype);
      % set(h, 'FaceVertexCData', 0.5*ones(length(distance),30));
      % set(h, 'FaceVertexCData', [0 0 1]);
      if strcmp(cfg.surftype, 'faces')
        set(h, 'edgealpha', 'interp')
        % set(h, 'edgealpha', 1)
        set(h, 'FaceVertexAlphaData', (1-distance).*(distance<0.1));
        % set(h, 'AlphaDataMapping', 'direct');
        set(h, 'facealpha', 'interp')
        % set(h, 'facealpha', 0)
        set(h, 'facecolor', 'none');
        % set(h, 'linestyle', 'none');
      end
    end
  end % plotspheres
  
  if  strcmp(cfg.plotbnd, 'yes') && ~isempty(vol) && isbem
    
    h = triplot(vol.bnd.pnt, vol.bnd.tri, [], cfg.surftype);
    % set(h, 'FaceVertexCData', 0.5*ones(length(distance),30));
    % set(h, 'FaceVertexCData', [0 0 1]);
    set(h, 'edgecolor', inner_skull_surface)
    % set(h, 'edgealpha', 1)
    if strcmp(cfg.surftype, 'faces')
      % set(h, 'AlphaDataMapping', 'direct');
      set(h, 'facealpha', 'interp')
      % set(h, 'facealpha', 0)
      set(h, 'facecolor', 'none');
      % set(h, 'linestyle', 'none');
    end
    
  end % plotbnd
  
  if strcmp(cfg.plotspherecenter, 'yes') && ~isempty(vol) && issphere
    plot3(vol.o(:,1), vol.o(:,2), vol.o(:,3), 'k.');
  end % plotspherecenter
  
  if strcmp(cfg.plotheadsurface, 'yes') && ~isempty(vol)
    % estimate the head surface from the spheres and gradiometers
    [pnt, tri] = headsurface(vol, sens);
    h = triplot(pnt, tri, [], cfg.surftype);
    set(h, 'edgecolor', cfg.surface_edgecolor);
    if strcmp(cfg.surftype, 'faces')
      set(h, 'facecolor', cfg.surface_facecolor);
      set(h, 'facealpha', cfg.surface_facealpha);
    end
  end % plotheadsurface
  
  if strcmp(cfg.plotlines, 'yes') && ismultisphere
    % first determine the indices of the relevant gradiometers.
    [sel_g, sel_v] = match_str(chan.label, vol.label);
    % create head-surface points from multisphere head-model.
    dir     = chan.pnt(sel_g,:) - vol.o(sel_v,:);
    dist    = sqrt(sum(dir.*dir,2));
    pnt0    = repmat((vol.r(sel_v)./dist),1,3).*dir + vol.o(sel_v,:);
    pnt1    = chan.pnt(sel_g,:);
    x = [pnt0(:,1) pnt1(:,1)]';
    y = [pnt0(:,2) pnt1(:,2)]';
    z = [pnt0(:,3) pnt1(:,3)]';
    line(x, y, z, 'color', 'm');
  end % plotlines
  
  if strcmp(cfg.plotfiducial, 'yes') && ~isempty(cfg.fiducial),
    fiduc = cfg.fiducial;
    plot3(fiduc(:,1), fiduc(:,2), fiduc(:,3), 'mo', 'lineWidth', 4);
  end
  
end % iseeg or ismeg

lighting gouraud
camlight
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
