function headmodel = ft_headmodel_duneuro(mesh, varargin)

% FT_HEADMODEL_DUNEURO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG. Different source models
% are implemented, including the St. Venant, the subtraction and partial
% integration model. This function takes as input a mesh with tetrahedral
% or hexahedral elements and corresponding conductivities and returns
% as output a volume conduction model which can be used to compute EEG
% leadfields.
%
% Use as
%   headmodel = ft_headmodel_duneuro(mesh,'conductivity', conductivities, ...)
%   headmodel = ft_headmodel_duneuro(mesh,'grid_filename', grid_filename, 'tensors_filename', tensors_filename, ...)
%
% Required input arguments should be specified in key-value pairs and have
% to include either
%   grid_filename   = string, filename for grid
%   tensors_filename= string, filename for conductivities
%
% or
%   conductivity    = vector, conductivity values for tissues
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   type            =  string, 'fitted' (default)
%   solver_type     =  string, 'cg' (default)
%   electrodes      =  string, 'closest_subentity_center' (default)
%   subentities     =  string, e.g. '1 2 3' (default) or '3'
%   forward         =  string, 'venant' (default), 'partial_integration' or 'subtraction'
%   intorderadd     =  string, e.g. '2' (default)
%   intorderadd_lb  =  string, e.g. '2' (default)
%   initialization  =  string, e.g. 'closest_vertex' (default)
%   numberOfMoments =  string, e.g. '3' (default)
%   referenceLength =  string, e.g. '20' (default)
%   relaxationFactor=  string, e.g. '1e-6' (default)
%   restrict        =  string, e.g. 'true' (default)
%   weightingExponent= string, e.g. '1' (default)
%   post_process    =  string, e.g. 'true' (default)
%   subtract_mean   =  string, e.g. 'true' (default)
%   reduction       =  string, e.g. '1e-10' (default)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD
%
% Copyright (C) 2017, Sophie Schrader
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
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.

%
% $Id$


% get the optional arguments and set defaults
type            = ft_getopt(varargin, 'type', 'fitted');
solver_type     = ft_getopt(varargin, 'solver_type', 'cg');
grid_filename   = ft_getopt(varargin, 'grid_filename');
tensors_filename= ft_getopt(varargin, 'tensors_filename');
conductivity    = ft_getopt(varargin, 'conductivity');
electrodes      =  ft_getopt(varargin, 'electrodes', 'closest_subentity_center');
subentities     =  ft_getopt(varargin, 'subentities', '1 2 3');
forward         =  ft_getopt(varargin, 'forward', 'venant');
intorderadd     =  ft_getopt(varargin, 'intorderadd', '2');
intorderadd_lb  =  ft_getopt(varargin, 'intorderadd_lb', '2');
initialization  =  ft_getopt(varargin, 'initialization', 'closest_vertex');
numberOfMoments =  ft_getopt(varargin, 'numberOfMoments', '3');
referenceLength =  ft_getopt(varargin, 'referenceLength', '20');
relaxationFactor=  ft_getopt(varargin, 'relaxationFactor', '1e-6');
restrict        =  ft_getopt(varargin, 'restrict', 'true');
weightingExponent =  ft_getopt(varargin, 'weightingExponent', '1');
post_process    =  ft_getopt(varargin, 'post_process', 'true');
subtract_mean   =  ft_getopt(varargin, 'subtract_mean', 'true');
reduction       =  ft_getopt(varargin, 'reduction', '1e-10');

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

if isfield(mesh,'tissue') && length(conductivity)<length(unique(headmodel.tissue))
  error('Wrong conductivity information!')
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
cfg = [];
cfg.type = type;
cfg.solver_type = solver_type;

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
  
  %make sure labels start at 0
  utissue = sort(unique(headmodel.tissue));
  for i = 1:size(headmodel.tissue,1)
    headmodel.tissue(i,1) =  find(utissue == headmodel.tissue(i,1));
  end
  cfg.volume_conductor.tensors.labels = uint64(headmodel.tissue -1);
  cfg.volume_conductor.tensors.conductivities = conductivity;
end

headmodel.driver = duneuro_meeg(cfg);

cfg                 = [];
cfg.type            = forward;

% optional arguments for subtraction and for Venant approach
cfg.intorderadd     = intorderadd;
cfg.intorderadd_lb  = intorderadd_lb;
cfg.initialization  = initialization;
cfg.numberOfMoments = numberOfMoments;
cfg.referenceLength = referenceLength;
cfg.relaxationFactor = relaxationFactor;
cfg.restrict        = restrict;
cfg.weightingExponent = weightingExponent;

%set source model
headmodel.driver.set_source_model(cfg);

headmodel.type = 'duneuro';
headmodel.electrodes = electrodes;
headmodel.subentities = subentities;
headmodel.post_process = post_process;
headmodel.subtract_mean = subtract_mean;
headmodel.reduction = reduction;

end
