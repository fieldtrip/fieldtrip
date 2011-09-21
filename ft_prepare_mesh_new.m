function bnd = ft_prepare_mesh(cfg, mri)

% FT_PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in world coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg, mri)
%
% Configuration options:
%   cfg.interactive     = 'no' (default) or 'yes' (manual interaction)
%   cfg.tissue          = list with segmentation values corresponding with each compartment
%   cfg.numvertices     = vector, length equal cfg.tissue.  e.g. [2000 1000 800];
%   cfg.downsample      = integer (1,2, ...) defines the level of refinement of the mri data
%   cfg.unit            = e.g. 'mm'
%   cfg.headshape       = (optional) a filename containing headshape, a Nx3 matrix with surface
%                         points, or a structure with a single or multiple boundaries
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% Example use:
%   mri=ft_read_mri('Subject01.mri');
%   cfg=[];
%   cfg.output = {'scalp', 'skull', 'brain'};
%   segment    = ft_volumesegment(cfg, mri);
%   scalp = (segment.scalp);
%   skull = 2*(segment.skull);
%   brain = 3*(segment.brain);
%   segment.seg=scalp+skull+brain;
%   cfg=[];
%   cfg.tissue      = [1 2 3];
%   cfg.numvertices = [2000 1000 800];
%   cfg.unit        = segment.unit;
%   bnd = ft_prepare_mesh(cfg, segment);

% Copyrights (C) 2009, Cristiano Micheli & Robert Oostenveld
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
% $Id: ft_prepare_mesh.m 4119 2011-09-06 14:05:57Z johzum $

cfg = ft_checkconfig(cfg, 'forbidden', 'numcompartments');

% set the defaults
downsample   = ft_getopt('downsample',cfg,1);
resolution   = ft_getopt('resolution',cfg,[]);
tissue       = ft_getopt('tissue',cfg,[]);
numvertices  = ft_getopt('numvertices',cfg,[]);
interactive  = ft_getopt('interactive',cfg,'no');
headshape    = ft_getopt('headshape',cfg,[]);
inputfile    = ft_getopt('inputfile',cfg,[]);
outputfile   = ft_getopt('outputfile',cfg,[]);
unit         = ft_getopt('unit',cfg,'mm');

if ~isempty(headshape) && isa(headshape, 'config')
  % convert the nested cmethodonfig-object back into a normal structure
  headshape = struct(headshape);
end

% load optional given inputfile like already segmented volume 
hasdata = (nargin>1);
if      hasdata && ~isempty(inputfile)
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
elseif  hasdata &&  isempty(inputfile)
  % this is ok
elseif ~hasdata && ~isempty(inputfile)
  % the input data should be read from file
  mri = loadvar(inputfile, 'mri');
elseif ~hasdata &&  isempty(inputfile)
  mri = [];
end

if isempty(headshape) && hasdata
  basedonseg        = isfield(mri, 'transform') && issegmentation(mri);
  basedonmri        = isfield(mri, 'transform') && ~basedonseg && isfield(mri, {'anatomy', 'dim'});
  basedonvol        = isfield(mri, 'bnd');
  basedonsphere     = isfield(mri,'r') && isfield(mri,'o');
  basedonheadshape  = 0;
  if ~(basedonseg+basedonmri+basedonvol+basedonsphere+basedonheadshape)
    error('unknown input type')
  end
elseif ~isempty(headshape) && ~hasdata
  basedonseg        = 0;
  basedonmri        = 0;
  basedonvol        = 0;
  basedonsphere     = 0;
  basedonheadshape  = 1;
elseif isempty(headshape) && ~hasdata
  error('no data available')  
else
  error('inconsistent configuration, cfg.headshape should not be used in combination with an mri input')
end

if basedonseg || basedonmri
  if downsample~=1
    % optionally downsample the anatomical MRI and/or the tissue segmentation
    fprintf('downsampling the volume by a factor of %d\n',downsample);
    cfg = [];
    cfg.outputfile = outputfile;
    cfg.downsample = downsample;
    mri = ft_volumedownsample(cfg, mri);
  end
end

if strcmp(interactive, 'yes')
  fprintf('using the manual approach\n');
  bnd = prepare_mesh_manual(cfg, mri);
  
elseif basedonmri
  fprintf('using the mri approach\n');
  mri = ft_datatype_volume(mri);
  % reslicing if necessary
  if ~isempty(resolution)
    fprintf('reslicing with homogeneous voxel resolution of %d %s\n',resolution,unit);
    cfg = [];
    cfg.resolution = resolution; 
    cfg.dim = mri.dim;
    mri_r   = ft_volumereslice(cfg, mri);
  end
  % segmenting the volume
  fprintf('segmenting into scalp/skull/brain compartments, this may take a while...\n');
  cfg = [];
  cfg.output = {'scalp','skull','brain'};
  [mri_s] = ft_volumesegment(cfg, mri_r);
  
  % A check/fix on the segmented volumes. It requires the image processing toolbox!
  fprintf('checking volumes\n');
  [mri_s.scalp] = fixhollow(mri_s.scalp);
  [mri_s.skull] = fixhollow(mri_s.skull);
  [mri_s.skull] = fixhollow(mri_s.skull);
  
  % preparing the meshes from the segmented compartments
  fprintf('preparing the meshes\n');
  tissuelabel = {'scalp','skull','brain'};
  bnd = prepare_mesh_segmentation_new(mri,'tissuelabel',tissuelabel);
  
elseif basedonseg 
  fprintf('using the segmentation approach\n');
  cfg = [];
  cfg.tissuelabel = cfg.segment;
  bnd = prepare_mesh_segmentation(cfg, mri);
  
elseif basedonheadshape
  fprintf('using the head shape to construct a triangulated mesh\n');
  cfg = [];
  cfg.headshape = headshape;
  bnd = prepare_mesh_headshape(cfg);
 
elseif basedonvol
  fprintf('using the mesh specified in the input volume conductor\n');
  bnd = mri.bnd;
  
elseif basedonsphere
  vol = mri;
  
  if isempty(numvertices)
    fprintf('using the mesh specified by icosaedron162\n');
    [pnt,tri] = icosahedron162;
  elseif any(numvertices==[42 162 642 2562])
    sprintf('using the mesh specified by icosaedron%d\n',numvertices);
    eval(['[pnt,tri] = icosahedron' num2str(numvertices) ';']);
  else
    [pnt, tri] = msphere(numvertices);
    sprintf('using the mesh specified by msphere with %d vertices\n',size(pnt,1));
  end
  
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
  end
  
else
  error('unsupported cfg.method and/or input')
end

% ensure that the vertices and triangles are double precision, otherwise the bemcp mex files will crash
for i=1:length(bnd)
  bnd(i).pnt = double(bnd(i).pnt);
  bnd(i).tri = double(bnd(i).tri);
end

% the output data should be saved to a MATLAB file
if ~isempty(outputfile)
  savevar(outputfile, 'data', bnd); % use the variable name "data" in the output file
end

function cmprtmnt = fixhollow(cmprtmnt)
% checks if the compartment is hollow
[~,N] = bwlabeln(cmprtmnt);
if N>2
% ...and tries to fix it
  cmprtmnt = imfill(cmprtmnt,'holes');
end

function res = issegmentation(mri,cfg)
res = false;
res = res || any(isfield(mri, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}));
% checks for existence of fields declared in the cfg.segment option
if isfield(cfg,'segment')
  for i=1:numel(cfg.segment)
    res = res || isfield(mri,cfg.segment{i});
  end
end
