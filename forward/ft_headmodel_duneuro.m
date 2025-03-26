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
%   meg_enablecache = logical, e.g. 'false (default)

ft_hastoolbox('duneuro', 1);

grid_filename    = ft_getopt(varargin, 'grid_filename');
tensors_filename = ft_getopt(varargin, 'tensors_filename');
conductivity     = ft_getopt(varargin, 'conductivity');
duneuro_settings = ft_getopt(varargin, 'duneuro_settings');

% get the optional arguments and set the defaults
duneuro_settings.type             = ft_getopt(duneuro_settings, 'type',             'fitted');
duneuro_settings.solver_type      = ft_getopt(duneuro_settings, 'solver_type',      'cg');
duneuro_settings.electrodes       = ft_getopt(duneuro_settings, 'electrodes',       'closest_subentity_center');
duneuro_settings.subentities      = ft_getopt(duneuro_settings, 'subentities',      '1 2 3');
duneuro_settings.forward          = ft_getopt(duneuro_settings, 'forward',          'venant');
duneuro_settings.intorderadd      = ft_getopt(duneuro_settings, 'intorderadd',      '2');
duneuro_settings.intorderadd_lb   = ft_getopt(duneuro_settings, 'intorderadd_lb',   '2');
duneuro_settings.initialization   = ft_getopt(duneuro_settings,' initialization',   'closest_vertex');
duneuro_settings.numberOfMoments  = ft_getopt(duneuro_settings, 'numberOfMoments',  '3');
duneuro_settings.referenceLength  = ft_getopt(duneuro_settings, 'referenceLength',  '20');
duneuro_settings.relaxationFactor = ft_getopt(duneuro_settings, 'relaxationFactor', '1e-6');
duneuro_settings.restrict         = ft_getopt(duneuro_settings, 'restrict',         'true');
duneuro_settings.weightingExponent = ft_getopt(duneuro_settings, 'weightingExponent', '1');
duneuro_settings.post_process      = ft_getopt(duneuro_settings, 'post_process',  'true');
duneuro_settings.subtract_mean     = ft_getopt(duneuro_settings, 'subtract_mean', 'true');
duneuro_settings.reduction         = ft_getopt(duneuro_settings, 'reduction', '1e-15');
duneuro_settings.intorderadd_meg   = ft_getopt(duneuro_settings, 'intorderadd_meg', '2');
duneuro_settings.mixedMoments      = ft_getopt(duneuro_settings, 'mixedMoments', 'true');
duneuro_settings.meg_type          = ft_getopt(duneuro_settings, 'meg_type', 'physical');
duneuro_settings.meg_enablecache   = ft_getopt(duneuro_settings, 'meg_enablecache', 'false');

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

% create cfg for the generation of the driver object (CentOS support only)
cfg                 = [];
cfg.type            = duneuro_settings.type;
cfg.solver_type     = duneuro_settings.solver_type;
cfg.meg.intorderadd = duneuro_settings.intorderadd_meg;
cfg.meg.type        = duneuro_settings.meg_type;
cfg.meg.cache.enable = duneuro_settings.meg_enablecache;

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

try
  % this only works (at least in my (=JM's) hands on a CentOS7 Linux machine)
  headmodel.driver = duneuro_meeg(cfg);
catch
  warning('To generate a valid Mex File, please compile duneuro-matlab on your machine following the instructions on the <a href = "https://gitlab.dune-project.org/duneuro/duneuro/-/wikis/Installation-instructions">DUNEuro Wiki page</a>.');
  return
end
headmodel.type = 'duneuro';
headmodel = copyfields(duneuro_settings, headmodel, {'forward' 'electrodes' 'subentities' 'intorderadd' 'intorderadd_lb' 'initialization', ...
              'numberOfMoments' 'referenceLength' 'relaxationFactor' 'restrict' 'weightingExponent' 'mixedMoments' 'post_process', ...
              'subtract_mean' 'reduction'});
            