function bnd = prepare_mesh(cfg, mri)

% PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in world coordinates.
%
% Use as
%   bnd = prepare_mesh(cfg, mri)
%
% Configuration options:
%   cfg.method          = 'segmentation' or 'manual'
%   cfg.tissue          = list with segmentation values corresponding with each compartment
%   cfg.downsample      = integer (1,2, ...) defines the level of refinement of the mri data
%   cfg.headshape     = a filename containing headshape, a Nx3 matrix with surface
%                       points, or a structure with a single or multiple boundaries
%
% Example use:
%   mri = read_mri('Subject01.mri');
%   cfg            = [];
%   cfg.method     = 'manual';
%   cfg.downsample = 2;
%   bnd = prepare_mesh(cfg, mri);

% Copyrights (C) 2009, Cristiano Micheli & Robert Oostenveld
%
% $Log: prepare_mesh.m,v $
% Revision 1.11  2009/07/29 06:43:16  roboos
% fixed basedonseg for headshape input (thanks to Vladimir)
%
% Revision 1.10  2009/07/16 09:00:51  crimic
% added choice of  multiple mesh methods and fixed a small typo
%
% Revision 1.9  2009/06/17 13:38:41  roboos
% cleaned up handling of method=manual
%
% Revision 1.8  2009/06/15 14:01:08  roboos
% minor updates in documentation and default cfg
%
% Revision 1.7  2009/06/03 12:07:11  crimic
% changed help
%
% Revision 1.6  2009/06/03 12:05:06  crimic
% added cfg.numcompartments option as input for automatic segmentation
%
% Revision 1.5  2009/06/02 10:18:39  crimic
% minor changes
%
% Revision 1.4  2009/05/14 19:23:33  roboos
% small cleanup, nothing functionally changed
%
% Revision 1.3  2009/05/06 16:09:17  roboos
% renamed gui_mesh into prepare_mesh_manual and moved to private
% some cleanup of the code
%
% Revision 1.2  2009/05/06 08:46:29  crimic
% First implementation
%

cfg = checkconfig(cfg, 'forbidden', 'numcompartments');

% set the defaults
if ~isfield(cfg, 'downsample'),      cfg.downsample = 1;         end
if ~isfield(cfg, 'tissue'),          cfg.tissue = [];            end
if ~isfield(cfg, 'numvertices'),     cfg.numvertices = [];       end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested cmethodonfig-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% there are three types of input possible
if nargin>1 && (~isfield(cfg,'headshape') || isempty(cfg.headshape))
  basedonseg        = isfield(mri, 'transform') && any(isfield(mri, {'seg', 'csf', 'white', 'gray'}));
  basedonmri        = isfield(mri, 'transform') && ~basedonseg;
  basedonvol        = isfield(mri, 'bnd');
  basedonheadshape  = 0;
elseif nargin==1 && isfield(cfg,'headshape') && ~isempty(cfg.headshape)
  basedonseg        = 0;
  basedonmri        = 0;
  basedonvol        = 0;
  basedonheadshape  = 1;
else
  error('inconsistent configuration, cfg.headshape should not be used in combination with an mri input')
end

if basedonseg || basedonmri
  % optionally downsample the anatomical MRI and/or the tissue segmentation
  tmpcfg = [];
  tmpcfg.downsample = cfg.downsample;
  mri = volumedownsample(tmpcfg, mri);
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
  
else
  error('unsupported cfg.method and/or input')
end

% ensure that the vertices and triangles are double precision, otherwise the bemcp mex files will crash
for i=1:length(bnd)
  bnd(i).pnt = double(bnd(i).pnt);
  bnd(i).tri = double(bnd(i).tri);
end

