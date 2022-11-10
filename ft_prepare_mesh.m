function [mesh, cfg] = ft_prepare_mesh(cfg, data)

% FT_PREPARE_MESH creates a triangulated surface mesh or tetrahedral/hexahedral
% volume mesh that can be used as geometrical description for a volume conduction
% model. The mesh can either be created manually from anatomical MRI data or can be
% generated starting from a segmented MRI. This function can also be used to create a
% cortex hull, i.e. the smoothed envelope around the pial surface created by
% freesurfer.
%
% Use as
%   mesh = ft_prepare_mesh(cfg)
%   mesh = ft_prepare_mesh(cfg, mri)
%   mesh = ft_prepare_mesh(cfg, seg)
% where the mri input argument is the result from FT_READ_MRI, FT_VOLUMEREALIGN or
% FT_VOLUMERESLICE and the seg input argument is from FT_VOLUMESEGMENT. If you
% specify an anatomical MRI, it will be segmented on the fly.
%
% The cfg argument is a structure that can contain:
%   cfg.method      = string, can be 'interactive', 'projectmesh', 'iso2mesh', 'isosurface',
%                     'headshape', 'hexahedral', 'tetrahedral', 'cortexhull' or 'fittemplate'
%   cfg.tissue      = cell-array with strings representing the tissue types, or numeric vector with integer values
%   cfg.numvertices = numeric vector, should have same number of elements as the number of tissues
%
% When providing an anatomical MRI or a segmentation, you should specify
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling), see FT_VOLUMEDOWNSAMPLE
%   cfg.spmversion  = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
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
%   mesh            = ft_prepare_mesh(cfg, segmentation);
%
%   cfg             = [];
%   cfg.method      = 'cortexhull';
%   cfg.headshape   = '/path/to/surf/lh.pial';
%   cfg.fshome      = '/path/to/freesurfer dir';
%   cortex_hull     = ft_prepare_mesh(cfg);
%
% See also FT_VOLUMESEGMENT, FT_PREPARE_HEADMODEL, FT_PLOT_MESH

% Undocumented functionality: at this moment it allows for either
%   mesh = ft_prepare_mesh(cfg)             or
%   mesh = ft_prepare_mesh(cfg, headmodel)
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
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% we cannot use nargin, because the data might have been loaded from cfg.inputfile
hasdata = exist('data', 'var');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'numcompartments', 'outputfile', 'sourceunits', 'mriunits'});

% get the options
cfg.downsample  = ft_getopt(cfg, 'downsample', 1);            % default is no downsampling
cfg.numvertices = ft_getopt(cfg, 'numvertices', 3000, true);  % set the default, [] is also a meaningful value
cfg.smooth      = ft_getopt(cfg, 'smooth');                   % no default
cfg.spmversion  = ft_getopt(cfg, 'spmversion', 'spm12');

% Translate the input options in the appropriate default for cfg.method
if isfield(cfg, 'headshape') && ~isempty(cfg.headshape)
  cfg.method = ft_getopt(cfg, 'method', 'headshape');
elseif hasdata && ~strcmp(ft_headmodeltype(data), 'unknown')
  cfg.method = ft_getopt(cfg, 'method', ft_headmodeltype(data));
elseif hasdata
  cfg.method = ft_getopt(cfg, 'method', 'projectmesh');
else
  cfg.method = ft_getopt(cfg, 'method', []);
end

if hasdata && cfg.downsample~=1
  % optionally downsample the anatomical volume and/or tissue segmentations
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  data = ft_volumedownsample(tmpcfg, data);
  % restore the provenance information and put back cfg.smooth
  tmpsmooth = cfg.smooth;
  [cfg, data] = rollback_provenance(cfg, data);
  cfg.smooth = tmpsmooth;
end

switch cfg.method
  case 'interactive'
    % this makes sense with a non-segmented MRI as input
    % call the corresponding helper function
    mesh = prepare_mesh_manual(cfg, data);
    
  case {'projectmesh', 'iso2mesh', 'isosurface'}
    % this makes sense with a segmented MRI as input
    % call the corresponding helper function
    mesh = prepare_mesh_segmentation(cfg, data);
    
  case 'headshape'
    % call the corresponding helper function
    mesh = prepare_mesh_headshape(cfg);
    
  case 'hexahedral'
    % the MRI is assumed to contain a segmentation
    % call the corresponding helper function
    mesh = prepare_mesh_hexahedral(cfg, data);
    
  case 'tetrahedral'
    % the input data is assumed to contain a segmentation
    % call the corresponding helper function
    mesh = prepare_mesh_tetrahedral(cfg, data);
    
  case {'singlesphere' 'concentricspheres'}
    % the input data is assumed to contain a spherical head model
    data = ft_datatype_headmodel(data);   % ensure that it is consistent and up-to-date
    data = ft_determine_units(data);      % ensure that it has units
    
    if length(data.r)>1 && numel(cfg.numvertices)==1
      % use the same number of vertices for each tissue
      cfg.numvertices = repmat(cfg.numvertices, length(data.r));
    end
    
    mesh = [];
    for i=1:length(data.r)
      ft_info('triangulating sphere %d in the volume conductor\n', i);
      [pos, tri] = mesh_sphere(cfg.numvertices(i));
      mesh(i).pos(:,1) = pos(:,1)*data.r(i) + data.o(1);
      mesh(i).pos(:,2) = pos(:,2)*data.r(i) + data.o(2);
      mesh(i).pos(:,3) = pos(:,3)*data.r(i) + data.o(3);
      mesh(i).tri = tri;
    end
    
  case 'cortexhull'
    mesh = prepare_mesh_cortexhull(cfg);
    
  case 'fittemplate'
    M    = prepare_mesh_fittemplate(cfg.headshape.pos, cfg.template.pos);
    mesh = ft_transform_geometry(M, data);
    
  otherwise
    ft_error('unsupported cfg.method')
end

% copy the geometrical units from the input to the output
if ~isfield(mesh, 'unit') && hasdata && isfield(data, 'unit')
  for i=1:numel(mesh)
    mesh(i).unit = data.unit;
  end
elseif ~isfield(mesh, 'unit')
  mesh = ft_determine_units(mesh);
end

% copy the coordinate system from the input to the output
if ~isfield(mesh, 'coordsys') && hasdata && isfield(data, 'coordsys')
  for i=1:numel(mesh)
    mesh(i).coordsys = data.coordsys;
  end
end

% smooth the mesh
if ~isempty(cfg.smooth)
  tmpcfg = keepfields(cfg, {'smooth', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  tmpcfg.numvertices = []; % the number of vertices should not be changed
  tmpcfg.headshape = mesh;
  mesh = prepare_mesh_headshape(tmpcfg);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance mesh
ft_postamble history    mesh
