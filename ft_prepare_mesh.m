function [bnd, cfg] = ft_prepare_mesh(cfg, mri)

% FT_PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. FT_PREPARE_MESH can be used
% to create a cortex hull, i.e. the smoothed envelope around the pial
% surface created by freesurfer. The result is a bnd structure which
% contains the information about all segmented surfaces related to mri
% sand are expressed in world coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg)
%   bnd = ft_prepare_mesh(cfg, mri)
%   bnd = ft_prepare_mesh(cfg, seg)
%
% Configuration options:
%   cfg.method      = string, can be 'interactive', 'projectmesh', 'iso2mesh', 'isosurface',
%                     'headshape', 'hexahedral', 'tetrahedral','cortexhull', 'fittemplate'
%   cfg.tissue      = cell-array with tissue types or numeric vector with integer values
%   cfg.numvertices = numeric vector, should have same number of elements as cfg.tissue
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling), see FT_VOLUMEDOWNSAMPLE
%   cfg.spmversion  = string, 'spm2', 'spm8', 'spm12' (default = 'spm8')
%
% For method 'headshape' you should specify
%   cfg.headshape   = a filename containing headshape, a Nx3 matrix with surface
%                     points, or a structure with a single or multiple boundaries
%
% For method 'cortexhull' you should not give input data, but specify
%   cfg.headshape   = string, filename containing the pial surface computed by freesurfer recon-all
%
% For method 'fittemplate' you should specify
%   cfg.headshape   = a filename containing headshape
%   cfg.template    = a filename containing headshape
% With this method you are fitting the headshape from the configuration to the template; 
% the resulting affine transformation is applied to the input mesh (or set of meshes), 
% which is subsequently returned as output variable.
%
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
%   cfg             = [];
%   cfg.method      = 'cortexhull';
%   cfg.headshape   = '/path/to/surf/lh.pial';
%   cfg.fshome      = '/path/to/freesurfer dir';
%   cortex_hull     = ft_prepare_mesh(cfg);
%
% See also FT_VOLUMESEGMENT, FT_PREPARE_HEADMODEL, FT_PLOT_MESH

% Undocumented functionality: at this moment it allows for either
%   bnd = ft_prepare_mesh(cfg)             or
%   bnd = ft_prepare_mesh(cfg, headmodel)
% but more consistent would be to specify a volume conduction model with
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%   cfg.headshape     = name of file containing the volume conduction model, see FT_READ_HEADMODEL
%
% Undocumented options, I have no clue why they exist
%   cfg.method = {'singlesphere' 'concentricspheres' 'localspheres'}

% Copyrights (C) 2009-2012, Robert Oostenveld & Cristiano Micheli
% Copyrights (C) 2012-2013, Robert Oostenveld & Lilla Magyari
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar mri
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% we cannot use nargin, because the data might have been loaded from cfg.inputfile
hasdata = exist('mri', 'var');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'numcompartments', 'outputfile', 'sourceunits', 'mriunits'});

% get the options
cfg.downsample  = ft_getopt(cfg, 'downsample', 1); % default is no downsampling
cfg.numvertices = ft_getopt(cfg, 'numvertices');   % no default
cfg.smooth      = ft_getopt(cfg, 'smooth');        % no default
cfg.spmversion  = ft_getopt(cfg, 'spmversion', 'spm8');

% Translate the input options in the appropriate default for cfg.method
if isfield(cfg, 'headshape') && ~isempty(cfg.headshape)
  cfg.method = ft_getopt(cfg, 'method', 'headshape');
elseif hasdata && ~strcmp(ft_headmodeltype(mri), 'unknown')
  cfg.method = ft_getopt(cfg, 'method', ft_headmodeltype(mri));
elseif hasdata
  cfg.method = ft_getopt(cfg, 'method', 'projectmesh');
else
  cfg.method = ft_getopt(cfg, 'method', []);
end

if hasdata && cfg.downsample~=1
  % optionally downsample the anatomical volume and/or tissue segmentations
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo'});
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information and put back cfg.smooth
  tmpsmooth = cfg.smooth;
  [cfg, mri] = rollback_provenance(cfg, mri);
  cfg.smooth = tmpsmooth;
end

switch cfg.method
  case 'interactive'
    % this makes sense with a non-segmented MRI as input
    % call the corresponding helper function
    bnd = prepare_mesh_manual(cfg, mri);

  case {'projectmesh', 'iso2mesh', 'isosurface'}
    % this makes sense with a segmented MRI as input
    % call the corresponding helper function
    bnd = prepare_mesh_segmentation(cfg, mri);

  case 'headshape'
    % call the corresponding helper function
    bnd = prepare_mesh_headshape(cfg);

  case 'hexahedral'
    % the MRI is assumed to contain a segmentation
    % call the corresponding helper function
    bnd = prepare_mesh_hexahedral(cfg, mri);

  case 'tetrahedral'
    % the MRI is assumed to contain a segmentation
    % call the corresponding helper function
    bnd = prepare_mesh_tetrahedral(cfg, mri);

  case {'singlesphere' 'concentricspheres'}
    headmodel = mri;
    headmodel = ft_datatype_headmodel(headmodel);   % ensure that it is consistent and up-to-date
    headmodel = ft_determine_units(headmodel);      % ensure that it has units
    bnd = [];
    [pos, tri] = mesh_sphere(cfg.numvertices);
    for i=1:length(headmodel.r)
      ft_info('triangulating sphere %d in the volume conductor\n', i);
      bnd(i).pos(:,1) = pos(:,1)*headmodel.r(i) + headmodel.o(1);
      bnd(i).pos(:,2) = pos(:,2)*headmodel.r(i) + headmodel.o(2);
      bnd(i).pos(:,3) = pos(:,3)*headmodel.r(i) + headmodel.o(3);
      bnd(i).tri = tri;
    end

  case 'cortexhull'
    bnd = prepare_mesh_cortexhull(cfg);
  
  case 'fittemplate'  
    M   = prepare_mesh_fittemplate(cfg.headshape.pos,cfg.template.pos);
    orig.mri = mri;
    orig = ft_transform_geometry(M,orig);
    bnd = orig.mri;
    
  otherwise
    ft_error('unsupported cfg.method')
end

% copy the geometrical units from the input to the output
if ~isfield(bnd, 'unit') && hasdata && isfield(mri, 'unit')
  for i=1:numel(bnd)
    bnd(i).unit = mri.unit;
  end
elseif ~isfield(bnd, 'unit')
  bnd = ft_determine_units(bnd);
end

% copy the coordinate system from the input to the output
if ~isfield(bnd, 'coordsys') && hasdata && isfield(mri, 'coordsys')
  for i=1:numel(bnd)
    bnd(i).coordsys = mri.coordsys;
  end
end

% smooth the mesh
if ~isempty(cfg.smooth)
  cfg.headshape = bnd;
  cfg.numvertices = [];
  bnd = prepare_mesh_headshape(cfg);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance bnd
ft_postamble history    bnd
