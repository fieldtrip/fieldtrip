function headmodel = ft_headmodel_duneuro(mesh, varargin)

% FT_HEADMODEL_DUNEURO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG and MEG. Different source models
% are implemented, including the St. Venant, the subtraction and partial
% integration model. This function takes as input a mesh with tetrahedral
% or hexahedral elements and corresponding conductivities and returns
% as output a volume conduction model which can be used to compute EEG/MEG
% leadfields.
%
% Use as
%   headmodel = ft_headmodel_duneuro(mesh,'conductivity', conductivities, ...)
%   headmodel = ft_headmodel_duneuro(mesh,'grid_filename', grid_filename, 'tensors_filename', tensors_filename, ...)
%
% Required input arguments should be specified in key-value pairs and have
% to include either
%   grid_filename   = string, filename for grid in "msh" fileformat (see here: https://gmsh.info/doc/texinfo/gmsh.html#File-formats)
%   tensors_filename= string, filename for conductivities, txt file with conductivity values
%
% or
%   conductivity    = vector, conductivity values for tissues
%
% Optional input arguments are passed with
%   duneuro_settings = (optional) struct, which can contain the following fields
%
%   type            = string, 'fitted' (default)
%   solver_type     = string, 'cg' (default)
%   electrodes      = string, 'closest_subentity_center' (default)
%   subentities     = string, e.g. '1 2 3' (default) or '3'
%   forward         = string, 'venant' (default), 'partial_integration'
%   intorderadd     = string, e.g. '2' (default)
%   intorderadd_lb  = string, e.g. '2' (default)
%   initialization  = string, e.g. 'closest_vertex' (default)
%   numberOfMoments = string, e.g. '3' (default)
%   referenceLength = string, e.g. '20' (default)
%   relaxationFactor= string, e.g. '1e-6' (default)
%   restrict        = string, e.g. 'true' (default)
%   weightingExponent= string, e.g. '1' (default)
%   post_process    = string, e.g. 'true' (default)
%   subtract_mean   = string, e.g. 'true' (default)
%   reduction       = string, e.g. '1e-10' (default)
%   intorderadd_meg = integer, e.g.'0' (default)
%   mixedMoments    = logical, e.g. 'true' (default)
%   meg_type        = string, e.g. 'physical' (default)
%   meg_eneablecache= logical, e.g. 'false (default)

ft_hastoolbox('duneuro', 1);

grid_filename   = ft_getopt(varargin, 'grid_filename');
tensors_filename= ft_getopt(varargin, 'tensors_filename');
conductivity    = ft_getopt(varargin, 'conductivity');
duneuro_settings = ft_getopt(varargin, 'duneuro_settings');

% get the optional arguments and set defaults
if(isfield(duneuro_settings, 'type')) type = duneuro_settings.type;
else type = 'fitted';
end

if(isfield(duneuro_settings, 'solver_type')) solver_type = duneuro_settings.solver_type;
else solver_type = 'cg';
end

if(isfield(duneuro_settings, 'electrodes')) electrodes = duneuro_settings.electrodes;
else electrodes = 'closest_subentity_center';
end

if(isfield(duneuro_settings, 'subentities')) subentities = duneuro_settings.subentities;
else subentities = '1 2 3';
end

if(isfield(duneuro_settings, 'forward')) forward = duneuro_settings.forward;
else forward = 'venant';
end

if(isfield(duneuro_settings, 'intorderadd')) intorderadd = duneuro_settings.intorderadd;
else intorderadd = '2';
end

if(isfield(duneuro_settings, 'intorderadd_lb')) intorderadd_lb = duneuro_settings.intorderadd_lb;
else intorderadd_lb = '2';
end

if(isfield(duneuro_settings, 'initialization')) initialization = duneuro_settings.initialization;
else initialization = 'closest_vertex';
end

if(isfield(duneuro_settings, 'numberOfMoments')) numberOfMoments = duneuro_settings.numberOfMoments;
else numberOfMoments = '3';
end

if(isfield(duneuro_settings, 'referenceLength')) referenceLength = duneuro_settings.referenceLength;
else referenceLength = '20';
end

if(isfield(duneuro_settings, 'relaxationFactor')) relaxationFactor = duneuro_settings.relaxationFactor;
else relaxationFactor = '1e-6';
end

if(isfield(duneuro_settings, 'restrict')) restrict = duneuro_settings.restrict;
else restrict = 'true';
end

if(isfield(duneuro_settings, 'weightingExponent')) weightingExponent = duneuro_settings.weightingExponent;
else weightingExponent = '1';
end

if(isfield(duneuro_settings, 'post_process')) post_process = duneuro_settings.post_process;
else post_process = 'true';
end

if(isfield(duneuro_settings, 'subtract_mean')) subtract_mean = duneuro_settings.subtract_mean;
else subtract_mean = 'true';
end

if(isfield(duneuro_settings, 'reduction')) reduction = duneuro_settings.reduction;
else reduction = '1e-15';
end


if(isfield(duneuro_settings, 'intorderadd_meg')) intorderadd_meg = duneuro_settings.intorderadd_meg;
else intorderadd_meg = '2';
end

if(isfield(duneuro_settings, 'mixedMoments')) mixedMoments = duneuro_settings.mixedMoments;
else mixedMoments = 'true';
end

if(isfield(duneuro_settings, 'meg_type')) reduction = duneuro_settings.meg_type;
else meg_type = 'physical';
end

if(isfield(duneuro_settings, 'meg_eneablecache')) meg_eneablecache = duneuro_settings.meg_eneablecache;
else meg_eneablecache = 'false';
end

% start with an empty volume conductor
mesh = ft_datatype_parcellation(mesh);
headmodel = [];
if isfield(mesh,'pos')
  headmodel.pos = mesh.pos;
else
  error('Vertex field is required!')
end

if isfield(mesh,'tet')
  headmodel.tet = mesh.tet;
elseif isfield(mesh,'hex')
  headmodel.hex = mesh.hex;
  %%assume fieldtrip ordering, reorder to conform to DUNE standard
  headmodel.hex(:,3) = mesh.hex(:,4);
  headmodel.hex(:,4) = mesh.hex(:,3);
  headmodel.hex(:,7) = mesh.hex(:,8);
  headmodel.hex(:,8) = mesh.hex(:,7);
else
  error('Element field with tetrahedral or hexahedral information is required!')
end

if isfield(mesh,'tissue')
  headmodel.tissue = mesh.tissue;
else
  error('No tissue labels declared!')
end

if (isempty(grid_filename) || isempty(tensors_filename)) && (isempty(conductivity))
  error('Either a filename for the grid and a filename for the conductivity or a conductivity array must be supplied')
end

if ~isempty(conductivity)
  if (size(conductivity,1) ~= 1 && size(conductivity,1) ~= 9)  || size(conductivity,2) ~= length(unique(headmodel.tissue))
    error('Conductivity has wrong dimension!')
  end
end

if ~isfield(mesh,'tissuelabel')
  numlabels = size(unique(mesh.tissue),1);
  headmodel.tissuelabel = {};
  ulabel = unique(labels);
  for i = 1:numlabels
    headmodel.tissuelabel{i} = num2str(ulabel(i));
  end
else
  headmodel.tissuelabel = mesh.tissuelabel;
end

% create driver object
cfg                 = [];
cfg.type            = type;
cfg.solver_type     = solver_type;
cfg.meg.intorderadd = intorderadd_meg;
cfg.meg.type        = meg_type;
cfg.meg.cache.enable = meg_eneablecache;


if isfield(mesh,'tet')
  cfg.element_type = 'tetrahedron';
elseif isfield(mesh,'hex')
  cfg.element_type = 'hexahedron';
end

if(~isempty(grid_filename) && (~isempty(tensors_filename)))
  cfg.volume_conductor.grid.filename = grid_filename;
  cfg.volume_conductor.tensors.filename = tensors_filename;
else
  if isfield(headmodel,'tet')
    cfg.volume_conductor.grid.elements = (uint64(headmodel.tet))' - 1;
  elseif isfield(headmodel,'hex')
    cfg.volume_conductor.grid.elements = (uint64(headmodel.hex))' - 1;
  end
  cfg.volume_conductor.grid.nodes = headmodel.pos';
  
  %make sure labels start at 0 %TODO: faster for 6C anisotropy
  utissue = sort(unique(headmodel.tissue));
  for i = 1:size(headmodel.tissue,1)
    headmodel.tissue(i,1) =  find(utissue == headmodel.tissue(i,1));
  end
  cfg.volume_conductor.tensors.labels = uint64(headmodel.tissue -1);
  if(size(conductivity,1) == 1)
    cfg.volume_conductor.tensors.conductivities = conductivity;
  elseif (size(conductivity,1) == 9)
    cfg.volume_conductor.tensors.tensors = conductivity;
  end
end

headmodel.driver = duneuro_meeg(cfg);

headmodel.type              = 'duneuro';
headmodel.forward           = forward;
headmodel.electrodes        = electrodes;
headmodel.subentities       = subentities;
headmodel.intorderadd       = intorderadd;
headmodel.intorderadd_lb    = intorderadd_lb;
headmodel.initialization    = initialization;
headmodel.numberOfMoments   = numberOfMoments;
headmodel.referenceLength   = referenceLength;
headmodel.relaxationFactor  = relaxationFactor;
headmodel.restrict          = restrict;
headmodel.weightingExponent = weightingExponent;
headmodel.mixedMoments      = mixedMoments;
headmodel.post_process      = post_process;
headmodel.subtract_mean     = subtract_mean;
headmodel.reduction         = reduction;

end
