function [bnd, cfg] = ft_prepare_mesh(cfg, mri)

% FT_PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in world coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg, mri)
%   bnd = ft_prepare_mesh(cfg, seg)
%
% Configuration options:
%   cfg.interactive     = 'no' (default) or 'yes' (manual interaction)
%   cfg.tissue          = cell-array with tissue types or numeric vector with integer values
%   cfg.numvertices     = numeric vector, should have same number of elements as cfg.tissue
%   cfg.downsample      = integer number (default = 1, i.e. no downsampling), see FT_VOLUMEDOWNSAMPLE
%   cfg.headshape       = (optional) a filename containing headshape, a Nx3 matrix with surface
%                         points, or a structure with a single or multiple boundaries
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% Example
%   mri             = ft_read_mri('Subject01.mri');
%
%   cfg             = [];
%   cfg.output      = {'scalp', 'skull', 'brain'};
%   segmentation    = ft_volumesegment(cfg, mri);
%
%   cfg             = [];
%   cfg.tissue      = {'scalp', 'skull', 'brain'};
%   cfg.numvertices = [800, 1600, 2400];
%   bnd             = ft_prepare_mesh(cfg, segmentation);
%
% See also FT_VOLUMESEGMENT, FT_PREPARE_HEADMODEL, FT_PLOT_MESH

% Undocumented functionality: at this moment it allows for either
%   bnd = ft_prepare_mesh(cfg)             or
%   bnd = ft_prepare_mesh(cfg, headmodel)
% but more consistent would be to specify a volume conduction model with
%   cfg.vol           = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%   cfg.headshape     = name of file containing the volume conduction model, see FT_READ_VOL
%
% Undocumented options
%   cfg.method = hexahedral, tetrahedral, isosurface

% Copyrights (C) 2009-2012, Robert Oostenveld & Cristiano Micheli
% Copyrights (C) 2012-2013, Robert Oostenveld & Lilla Magyari
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

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar mri

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'numcompartments', 'outputfile', 'sourceunits', 'mriunits'});

% get the defaults
cfg.headshape    = ft_getopt(cfg, 'headshape');         % input option
cfg.interactive  = ft_getopt(cfg, 'interactive', 'no'); % to interact with the volume
% cfg.tissue       = ft_getopt(cfg, 'tissue');            % to perform the meshing on a specific tissue
cfg.numvertices  = ft_getopt(cfg, 'numvertices');       % resolution of the mesh
cfg.downsample   = ft_getopt(cfg, 'downsample', 1);

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% here we cannot use nargin, because the data might have been loaded from cfg.inputfile
hasdata = exist('mri', 'var');

if ~hasdata
  mri = [];
elseif ~ft_voltype(mri, 'unknown')
  % The input appears to be a headmodel. This is deprecated, but at this
  % moment (2012-09-28) we decided not to break the old functionality yet.
else
  mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'});
end

if hasdata
  % try to estimate the units, these will also be assigned to the output meshes
  mri = ft_convert_units(mri);
end

if hasdata
  % determine the type of input data
  basedonmri        = ft_datatype(mri, 'volume');
  basedonseg        = ft_datatype(mri, 'segmentation');
  basedonheadshape  = 0;
  basedonbnd        = isfield(mri, 'bnd');
  basedonsphere     = all(isfield(mri, {'r', 'o'}));
elseif isfield(cfg,'headshape') && ~isempty(cfg.headshape)
  % in absence of input data
  basedonmri       = false;
  basedonseg       = false;
  basedonheadshape = true;
  basedonbnd       = false;
  basedonsphere    = false;
else
  error('inconsistent configuration and input data');
end

if isfield(cfg, 'method') && strcmp(cfg.method, 'hexahedral')
  % the MRI is assumed to contain a segmentation, call the helper function
  bnd = prepare_mesh_hexahedral(cfg, mri); %should go fieldtrip/private
  % ensure that non of the other options gets executed
  basedonmri       = false;
  basedonseg       = false;
  basedonheadshape = false;
  basedonbnd       = false;
  basedonsphere    = false;
  
elseif isfield(cfg, 'method') && strcmp(cfg.method, 'tetrahedral')
  % the MRI is assumed to contain a segmentation, call the helper function
  bnd = prepare_mesh_tetrahedral(cfg, mri);
  % ensure that non of the other options gets executed
  basedonmri       = false;
  basedonseg       = false;
  basedonheadshape = false;
  basedonbnd       = false;
  basedonsphere    = false;
  
elseif basedonseg || basedonmri
  if all(isfield(mri, {'gray', 'white', 'csf'}))
    cfg.tissue      = ft_getopt(cfg, 'tissue', 'brain');
    cfg.numvertices = ft_getopt(cfg, 'numvertices', 3000);
  else
    cfg.tissue      = ft_getopt(cfg, 'tissue');
  end
  cfg = ft_checkconfig(cfg, 'required', {'tissue', 'numvertices'});
end

if (basedonseg || basedonmri) && cfg.downsample~=1
  % optionally downsample the anatomical MRI and/or the tissue segmentation
  tmpcfg = [];
  tmpcfg.downsample = cfg.downsample;
  mri = ft_volumedownsample(tmpcfg, mri);
end

if (basedonmri || basedonseg) && istrue(cfg.interactive)
  % this only makes sense with a (segmented) MRI as input
  fprintf('using the manual approach\n');
  bnd = prepare_mesh_manual(cfg, mri);
  
elseif basedonseg
  % FIXME this should be renamed to prepare_mesh_triangulation
  fprintf('using the segmentation approach\n');
  bnd = prepare_mesh_segmentation(cfg, mri);
  
elseif basedonmri && iscell(cfg.tissue) && all(isfield(mri, cfg.tissue))
  % the input is not detected as segmentation, but it does have all fields
  % that the user requests to have triangulated, so assume that it is a
  % segmentation after all
  fprintf('using the segmentation approach\n');
  bnd = prepare_mesh_segmentation(cfg, mri);
  
elseif basedonmri
  error('Unsegmented MRI only allowed in combination with cfg.interactive=yes')
  
elseif basedonheadshape
  fprintf('using the head shape to construct a triangulated mesh\n');
  bnd = prepare_mesh_headshape(cfg);
  
elseif basedonbnd
  fprintf('using the mesh specified in the input volume conductor\n');
  bnd = mri.bnd;
  
elseif basedonsphere
  fprintf('triangulating the sphere in the volume conductor\n');
  vol = mri;
  
  [pnt, tri] = makesphere(cfg.numvertices);
  
  
  switch ft_voltype(vol)
    case {'singlesphere' 'concentricspheres'}
      bnd = [];
      for i=1:length(vol.r)
        bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(1);
        bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(2);
        bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(3);
        bnd(i).tri = tri;
      end
    case 'localspheres'
      % FIXME this should be replaced by an outline of the head, see private/headsurface
      bnd = [];
      for i=1:length(vol.label)
        bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(i,1);
        bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(i,2);
        bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(i,3);
        bnd(i).tri = tri;
      end
  end
  
elseif isfield(cfg, 'method') && strcmp(cfg.method, 'hexahedral')
  % do nothing
elseif isfield(cfg, 'method') && strcmp(cfg.method, 'tetrahedral')
  % do nothing
else
  error('unsupported cfg.method and/or input')
end

% copy the geometrical units from the input to the output
if ~isfield(bnd, 'unit') && hasdata
  for i=1:numel(bnd)
    bnd(i).unit = mri.unit;
  end
elseif ~isfield(bnd, 'unit')
  bnd = ft_convert_units(bnd);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous mri
ft_postamble history bnd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pnt, tri] = makesphere(numvertices)

if isempty(numvertices)
  [pnt,tri] = icosahedron162;
  fprintf('using the mesh specified by icosaedron162\n');
elseif numvertices==42
  [pnt,tri] = icosahedron42;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==162
  [pnt,tri] = icosahedron162;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==642
  [pnt,tri] = icosahedron642;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==2562
  [pnt,tri] = icosahedron2562;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
else
  [pnt, tri] = msphere(numvertices);
  fprintf('using the mesh specified by msphere with %d vertices\n',size(pnt,1));
end
