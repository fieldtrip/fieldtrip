function [bnd, cfg] = ft_prepare_mesh_new(cfg, data)

% FT_PREPARE_MESH_NEW creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in head coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg, data)
%
% Configuration options:
%   cfg.interactive     = 'no' (default) or 'yes' (manual interaction)
%   cfg.tissue          = list with segmentation values/names corresponding with each compartment
%   cfg.numvertices     = vector, length equal cfg.tissue.  e.g. [2000 1000 800];
%   cfg.downsample      = (optional) integer (1,2, ...) defines the level of refinement of the mri data
%   cfg.headshape       = a filename containing headshape, a Nx3 matrix with surface
%                         points, or a structure with a single or multiple boundaries
%   cfg.unit            = e.g. 'mm'
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% Example use:
%
%   mri = ft_read_mri('Subject01.mri');
%   cfg = [];
%   cfg.output = {'scalp', 'skull', 'brain'};
%   cfg.numvertices = [2000 1000 800];
%   bnd = ft_prepare_mesh(cfg, segment);
%
% See also FT_PREPARE_CONCENTRICSPHERES, FT_PREPARE_LOCALSPHERES,
% FT_PREPARE_SINGLESHELL, FT_PREPARE_LEADFIELD, FT_PREPARE_MESH,
% FT_PREPARE_BEMMODEL

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
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', 'numcompartments', ...
                                       'method');

% mri defaults
resolution   = ft_getopt(cfg, 'resolution');  % for mri reslice
resdim       = ft_getopt(cfg, 'resdim');      % for mri reslice
downsample   = ft_getopt(cfg, 'downsample');  % for mri downsample
output       = ft_getopt(cfg, 'output');      % for mri segment
smooth       = ft_getopt(cfg, 'smooth');      % for mri segment
threshold    = ft_getopt(cfg, 'threshold');   % for mri segment

% ft_prepare_mesh specific options
headshape    = ft_getopt(cfg, 'headshape');         % input option
inputfile    = ft_getopt(cfg, 'inputfile');         % input option
outputfile   = ft_getopt(cfg, 'outputfile');        % output option
interactive  = ft_getopt(cfg, 'interactive', 'no'); % to interact with the volume
smoothseg    = ft_getopt(cfg, 'smoothseg', 0);      % to smooth the seg compartment
thresholdseg = ft_getopt(cfg, 'thresholdseg', 0);   % to threshold the seg compartment
tissue       = ft_getopt(cfg, 'tissue');            % to perform the meshing on a specific tissue
numvertices  = ft_getopt(cfg, 'numvertices');       % resolution of the mesh
interactive  = istrue(interactive);

% mesh points' unit options
transform    = ft_getopt(cfg, 'transform', eye(4));
unit         = ft_getopt(cfg, 'unit');              % target units

hasdata = exist('data', 'var');

if hasdata
  % if the target units are not explicitly given, try to estimate the data
  % units and use them for the meshes
  if ft_datatype(data,'volume')
    if ~isfield(data, 'unit'), data = ft_convert_units(data); end
  end
  if isempty(unit)
    if ft_datatype(data,'volume')
      mriunits    = data.unit;
      sourceunits = mriunits;
    else
      sourceunits = 'mm';
      mriunits    = 'mm';
    end
  end
else
  % this is a trick not to reascale in case of boundary input
  sourceunits   = 'mm';
  mriunits      = 'mm';
end

if isempty(headshape) && ~hasdata && ~hasinputfile
  error('no input data available')
end

if ~isempty(headshape) && isa(headshape, 'config')
  % convert the nested config-object back into a normal structure
  headshape = struct(headshape);
end

if hasdata
  basedonmri    = ft_datatype(data,'volume') && ~issegmentation(data,cfg);
  basedonrawseg = (~isstruct(data) && numel(size(data))==3);
  basedonmriseg = ~basedonrawseg && issegmentation(data,cfg);
  basedonvol    = isfield(data, 'bnd');
  basedonsphere = isfield(data,'r') && isfield(data,'o');
  basedonheadshape = 0;
else
  % in absence of data input
  basedonmri    = 0;
  basedonrawseg = 0;
  basedonmriseg = 0;
  basedonvol    = 0;
  basedonsphere = 0;
  basedonheadshape = ~isempty(headshape);
end

if ~interactive
  if (basedonmri+basedonrawseg+basedonmriseg+basedonvol+basedonsphere+basedonheadshape)>1
    error('inconsistent configuration, input data is ambiguous')
  end
else
  if hasdata || hasinputfile
    bnd = prepare_mesh_manual(cfg, data);
    return
  else
    error('you need a data structure')
  end
end

if basedonmri
  fprintf('using the mri approach\n');
  mri = ft_datatype(data,'volume');
  if isempty(resdim)
    resdim = mri.dim;
  end
  
  if ~isfield(mri, 'unit'), mri = ft_convert_units(mri); end
  
  % reslicing if necessary
  if ~isempty(resolution)
    fprintf('reslicing with homogeneous voxel resolution of %d %s\n',resolution,unit);
    cfg = [];
    cfg.resolution = resolution; % FIXME: what happens if this is in cm
    cfg.dim        = resdim;
    mri = ft_volumereslice(cfg, mri);
  end
  
  % segmenting the volume
  cfg = [];
  cfg.smooth     = smooth;
  cfg.threshold  = threshold;
  cfg.output     = output;
  cfg.downsample = downsample;
  mri = ft_volumesegment(cfg,mri);
  
  cfg = [];
  cfg.tissue = output;
  bnd = ft_prepare_mesh_new(cfg, mri);
  
elseif basedonmriseg
  if isempty(cfg.tissue) 
    [dum,tissue] = issegmentation(data,cfg);
  elseif isnumeric(cfg.tissue)
    error('tissue type should be a string')
  end

  fprintf('using the mri segmentation approach\n');
  cfg = [];
  cfg.tissue       = tissue;
  cfg.thresholdseg = thresholdseg;
  cfg.smoothseg    = smoothseg;
  cfg.numvertices  = numvertices;
  bnd = prepare_mesh_segmentation_new(cfg, data);
  
elseif basedonrawseg
  if ~isnumeric(cfg.tissue)
    error('tissue type should be a number')
  end
  fprintf('using the raw segmentation approach\n');
  cfg = [];
  cfg.tissue       = tissue;
  cfg.thresholdseg = thresholdseg;
  cfg.smoothseg    = smoothseg;
  cfg.numvertices  = numvertices;
  bnd = prepare_mesh_segmentation_new(cfg, data);
  
elseif basedonheadshape
  fprintf('using the head shape to construct a triangulated mesh\n');
  cfg = [];
  cfg.headshape   = headshape;
  cfg.numvertices = numvertices;
  bnd = prepare_mesh_headshape(cfg);
  
elseif basedonvol
  vol = data;
  fprintf('using the mesh specified in the input volume conductor\n');
  for i=1:numel(vol)
    tmpbnd(i) = vol(i).bnd;
  end
  cfg = [];
  cfg.headshape   = tmpbnd;
  cfg.numvertices = numvertices;
  bnd = prepare_mesh_headshape(cfg);
  
elseif basedonsphere
  vol = data;
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
        bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(i,1);
        bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(i,2);
        bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(i,3);
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
  
end

% voxel to head coordinates transformation
for i=1:length(bnd)
  % ensure that the vertices and triangles are double precision, otherwise
  % the bemcp mex files will crash
  bnd(i).pnt = double(bnd(i).pnt);
  bnd(i).tri = double(bnd(i).tri);
  % apply the coordinate transformation from voxel to head coordinates
  bnd(i).pnt = warp_apply(transform,bnd(i).pnt);
  % rescale
  bnd(i).pnt = scaleunit(sourceunits,mriunits,bnd(i).pnt);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res,fnames] = issegmentation(mri,cfg)
res = false;
% res = res || any(isfield(mri, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}));
names = {'segmentation', 'segment', 'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'};
fnames = {};
tmp = fieldnames(mri);

cnt = 1;
for i=1:numel(tmp)
  if ismember(tmp{i},names)
    fnames{cnt} = tmp{i};
    res = true; cnt = cnt +1;
  end
end

% checks for existence of fields declared in the cfg.tissue option
if isfield(cfg,'tissue')
  if ~isempty(cfg.tissue) && ~isnumeric(cfg.tissue)
    res = true;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pnt = scaleunit(sourceunits, mriunits, pnt)
% convert the MRI surface points into the same units as the source/gradiometer
scale = 1;
switch sourceunits
  case 'mm'
    scale = scale * 1000;
  case 'cm'
    scale = scale * 100;
  case 'dm'
    scale = scale * 10;
  case 'm'
    scale = scale * 1;
  otherwise
    error('unknown physical dimension in input ''unit''');
end
switch mriunits
  case 'mm'
    scale = scale / 1000;
  case 'cm'
    scale = scale / 100;
  case 'dm'
    scale = scale / 10;
  case 'm'
    scale = scale / 1;
  otherwise
    error('unknown physical dimension in mri.unit');
end
if scale~=1
  fprintf('converting MRI surface points from %s into %s\n', sourceunits, mriunits);
  pnt = pnt*(scale);
end
