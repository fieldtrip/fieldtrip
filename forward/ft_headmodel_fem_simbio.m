function vol = ft_headmodel_fem_simbio(seg,varargin)
% FT_HEADMODEL_FEM_SIMBIO reads a volume conduction model from a Vista .v
% file
%
%
% Use as
%   vol = ft_headmodel_fem_simbio(seg)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD, TEST_HEADMODEL_SIMBIO

% Copyright (C) 2011, Cristiano Micheli

wfmethod     = ft_getopt(varargin, 'wfmethod', 'cubes'); % wireframe method (default: cubes)
transform    = ft_getopt(varargin, 'transform', []);  % contains the tr. matrix from voxels to head coordinates
unit         = ft_getopt(varargin, 'unit', 'cm');     % contains the units of the transformation matrix (mm, cm, ...)
tissue       = ft_getopt(varargin, 'tissue', []);     % contains the labels of the tissues
tissueval    = ft_getopt(varargin, 'tissueval', []);  % contains the tissue values (an integer for each compartment)
tissuecond   = ft_getopt(varargin, 'tissuecond', []); % contains the tissue conductivities
tissueweight = ft_getopt(varargin, 'tissueweight', ones(1,numel(tissueval))); % contains the weigths for vgrid
% (how dense the wireframe of each tissue conpartment should be done)

% start with an empty structure
labels = []; % used to write the .par file (integers like 101,102,...)

% define a wireframe for each compartment (FEM grid)
wf = [];
if strcmp(wfmethod,'cubes')
  
  % write the segmented volume in a Vista format .v file
  [~,tname] = fileparts(tempname);
  MRfile = [tname '.v'];
  % write a vista volumetric file
  ft_write_volume(MRfile, seg,'dataformat', 'vista'); % calls write_vista_vol(size(seg), seg, MRfile);

  % write the materials file (assign tissue values to the elements of the FEM grid)
  % see tutorial: http://www.rheinahrcampus.de/~medsim/vgrid/manual.html
  [~,tname] = fileparts(tempname);
  materialsfile = [tname '.mat'];
  write_materials(materialsfile,tissueval,tissue,tissueweight);
  
  % Use vgrid to get the wireframe
  [~,tname] = fileparts(tempname);
  meshfile  = [tname '.v'];
  [~,tname] = fileparts(tempname);
  shfile    = [tname '.sh'];
  %FIXME:  use hastoolbox, change this in the future
  vroot = '/home/coherence/crimic/test/SimBio/fromJohannes/vgrid1.3.1/program/';
  efid  = fopen(shfile, 'w');
  fprintf(efid,'#!/usr/bin/env bash\n');
  % works with voxels
  fprintf(efid,[vroot 'vgrid -in ' MRfile ' -out ' meshfile ' -min 1 -max 1 -elem cube', ...
    ' -material ' materialsfile ' -smooth shift -shift 0.49 2>&1 > /dev/null\n']);
  fclose(efid);
  % execute command
  dos(sprintf('chmod +x %s', shfile));
  disp('Vgrid is writing the mesh file, this may take some time ...')
  stopwatch = tic;
  dos(['./' shfile]);
  disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
  
  % read the mesh points
  [shape] = ft_read_headshape(meshfile,'format','vista','unit',unit); % calls [nodes,elements,labels] = read_vista_mesh(meshfile);
  nodes    = shape.nd; % in voxel coordinates
  elements = shape.el;
  labels   = shape.labels;
  
  % assign wireframe (or FEM grid, or 3d mesh)
  wf.nd = nodes; 
  wf.el = elements;
  wf.labels = labels;
  
elseif strcmp(wfmethod,'tetra')
  % not yet impementes
else
  error('Unsupported method')
end

vol = [];
vol.wf        = wf;
vol.cond      = tissuecond;
vol.labels    = labels;
vol.transform = transform;
vol.unit      = unit;
vol.type      = 'simbio';
