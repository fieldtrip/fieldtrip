function headmodel = ft_headmodel_duneuro(mesh, varargin)

% FT_HEADMODEL_DUNEURO creates a volume conduction model of the head using the finite element method (FEM) for EEG and MEG.
% Different source models are implemented, including the St. Venant, the subtraction and partial integration model. This 
% function takes as input a mesh with tetrahedral or hexahedral elements and corresponding conductivities and returns
% as output a volume conduction model which can be used to compute EEG/MEG leadfields.
%
% Use as
%   headmodel = ft_headmodel_duneuro(mesh,'conductivity', conductivities, ...)
%   headmodel = ft_headmodel_duneuro(mesh,'grid_filename', grid_filename, 'tensors_filename', tensors_filename, ...)
%
% Required input arguments should be specified in key-value pairs and have to include either
%   grid_filename   = string, filename for grid in "msh" fileformat (see here: https://gmsh.info/doc/texinfo/gmsh.html#File-formats)
%   tensors_filename= string, filename for conductivities, txt file with conductivity values
% or
%   conductivity    = vector, conductivity values for tissues
% 
% if a pair of filenames is provided, the input mesh is not considered, but will be generated from the grid_filename
% 
% In addition, an optional struct with configuration options can be provided, which can specify the options
% related to the functional behavior of the duneuro software. See DUNEURO_DEFAULTS for the configureable options, and 
% their default values.

ft_hastoolbox('duneuro', 1);

grid_filename    = ft_getopt(varargin, 'grid_filename');
tensors_filename = ft_getopt(varargin, 'tensors_filename');
conductivity     = ft_getopt(varargin, 'conductivity');
duneuro_settings = ft_getopt(varargin, 'duneuro_settings');

duneuro_settings = duneuro_defaults(duneuro_settings); 

% start with an empty volume conductor
headmodel = [];

if (isempty(grid_filename) || isempty(tensors_filename)) && (isempty(conductivity))
  error('Either a filename for the grid and a filename for the conductivity or a conductivity array must be supplied')
end

% either the grid and tensor filenames need to be defined, or the mesh should contain the geometry + a non-empty conductivity 
% variable, usemesh has priority. if usefiles is used, element_type should be defined in the input parameters 
usemesh  = ~isempty(mesh);
usefiles = ~isempty(grid_filename) && ~isempty(tensors_filename);
assert(usemesh || usefiles);

if usemesh
  mesh = ft_datatype_parcellation(mesh);
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
end

% create cfg for the generation of the driver object
cfg                 = [];
cfg.type            = duneuro_settings.type;
cfg.solver_type     = duneuro_settings.solver_type;
cfg.meg.intorderadd = num2str(duneuro_settings.meg.intorderadd, '%d');
cfg.meg.type        = duneuro_settings.meg.type;
cfg.meg.cache.enable = bool2str(duneuro_settings.meg.enablecache);

if usemesh
  if isfield(mesh,'tet')
    cfg.element_type = 'tetrahedron';
  elseif isfield(mesh,'hex')
    cfg.element_type = 'hexahedron';
  end
elseif usefiles
  assert(~isempty(element_type));
  cfg.element_type = element_type;
end

if usefiles
  cfg.volume_conductor.grid.filename = grid_filename;
  cfg.volume_conductor.tensors.filename = tensors_filename;
elseif usemesh
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
  if size(conductivity,1) == 1
    cfg.volume_conductor.tensors.conductivities = conductivity;
  elseif size(conductivity,1) == 9
    cfg.volume_conductor.tensors.tensors = conductivity;
    cfg.useTensor = false;
  end
end

if ~duneuro_settings.bstflag
  try
    headmodel.driver = duneuro_meeg(cfg);
  catch
    warning('To generate a valid Mex File, please compile duneuro-matlab on your machine following the instructions on the <a href = "https://gitlab.dune-project.org/duneuro/duneuro/-/wikis/Installation-instructions">DUNEuro Wiki page</a>.');
    return
  end
else
  % prepare for the creation of a ini-file and a single shot computation,
  % this requires to write a grid_filename file and a tensors_filename
  % file, as well as a file that holds the sensors/electrodes, and a file
  % that holds the source model. The latter 2 are not available yet.
  filename = fullfile(duneuro_settings.outputpath, 'headmodel');
  duneuro_settings.grid_filename = duneuro_write_headmodel(filename, headmodel, duneuro_settings);

  
  if size(conductivity,1)==1
    filename = fullfile(duneuro_settings.outputpath, 'cond');
    duneuro_settings.tensors_filename = duneuro_write_conductivity(filename, 'conductivity', conductivity);
  end

end
headmodel.type    = 'duneuro';
headmodel.duneuro = duneuro_settings;
headmodel.duneuro.element_type = cfg.element_type;

function output = bool2str(input)

if input
  output = 'true';
else
  output = 'false';
end