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
%   cfg.method          = 'segmentation' or 'manual'
%   cfg.tissue          = list with segmentation values corresponding with each compartment
%   cfg.downsample      = integer (1,2, ...) defines the level of refinement of the mri data
%   cfg.headshape     = a filename containing headshape, a Nx3 matrix with surface
%                       points, or a structure with a single or multiple boundaries
%
% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk
%
% Example use:
%   mri            = ft_read_mri('Subject01.mri');
%   cfg            = [];
%   cfg.method     = 'manual';
%   cfg.downsample = 2;
%   bnd = prepare_mesh(cfg, mri);

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

cfg = checkconfig(cfg, 'forbidden', 'numcompartments');

% set the defaults
if ~isfield(cfg, 'downsample'),      cfg.downsample = 1;         end
if ~isfield(cfg, 'tissue'),          cfg.tissue = [];            end
if ~isfield(cfg, 'numvertices'),     cfg.numvertices = [];       end
if ~isfield(cfg, 'inputfile'),       cfg.inputfile = [];         end
if ~isfield(cfg, 'outputfile'),      cfg.outputfile = [];        end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested cmethodonfig-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% load optional given inputfile like already segmented volume 
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    mri = loadvar(cfg.inputfile, 'data');
  end
else
  mri = [];
end

if ~isfield(cfg,'headshape') || isempty(cfg.headshape)
  basedonseg        = isfield(mri, 'transform') && any(isfield(mri, {'seg', 'csf', 'white', 'gray'}));
  basedonmri        = isfield(mri, 'transform') && ~basedonseg;
  basedonvol        = isfield(mri, 'bnd');
  basedonsphere     = isfield(mri,'r') && isfield(mri,'o');
  basedonheadshape  = 0;
elseif isfield(cfg,'headshape') && ~isempty(cfg.headshape)
  basedonseg        = 0;
  basedonmri        = 0;
  basedonvol        = 0;
  basedonsphere     = 0;
  basedonheadshape  = 1;
else
  error('inconsistent configuration, cfg.headshape should not be used in combination with an mri input')
end

if basedonseg || basedonmri
  % optionally downsample the anatomical MRI and/or the tissue segmentation
  tmpcfg = [];
  tmpcfg.outputfile = cfg.outputfile;
  tmpcfg.downsample = cfg.downsample;
  mri = ft_volumedownsample(tmpcfg, mri);
end

if basedonseg
  fprintf('using the segmentation approach\n');
  bnd = prepare_mesh_segmentation(cfg, mri);
  
elseif basedonmri
  fprintf('using the manual approach\n');
  bnd = prepare_mesh_manual(cfg, mri);
  
elseif basedonheadshape
  fprintf('using the head shape to construct a triangulated mesh\n');
  bnd = prepare_mesh_headshape(cfg);
  
elseif basedonvol
  fprintf('using the mesh specified in the input volume conductor\n');
  bnd = mri.bnd;
  
elseif basedonsphere
  vol = mri;
  
  if isempty(cfg.numvertices)
    fprintf('using the mesh specified by icosaedron162\n');
    [pnt,tri] = icosahedron162;
  elseif any(cfg.numvertices==[42 162 642 2562])
    sprintf('using the mesh specified by icosaedron%d\n',cfg.numvertices);
    eval(['[pnt,tri] = icosahedron' num2str(cfg.numvertices) ';']);
  else
    [pnt, tri] = msphere(cfg.numvertices);
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
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', bnd); % use the variable name "data" in the output file
end
