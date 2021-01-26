function headmodel = ft_headmodel_simbio(mesh, varargin)

% FT_HEADMODEL_SIMBIO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG. This function takes
% as input a volumetric mesh (hexahedral or tetrahedral) and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This implements
%       ...
%
% Use as
%   headmodel = ft_headmodel_simbio(mesh,'conductivity', conductivities, ...)
%
% The mesh is given as a volumetric mesh, using ft_datatype_parcellation
%   mesh.pos = vertex positions
%   mesh.tet/mesh.hex = list of volume elements
%   mesh.tissue = tissue assignment for elements
%   mesh.tissuelabel = labels correspondig to tissues
%
% Required input arguments should be specified in key-value pairs and have
% to include
%   conductivity   = vector containing tissue conductivities using ordered
%                    corresponding to mesh.tissuelabel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this on Windows the following packages are necessary:
%
% Microsoft Visual C++ 2008 Redistributable
%
% Intel Visual Fortran Redistributables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% $Id$

ft_hastoolbox('simbio', 1);

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');

% start with an empty volume conductor
mesh = ft_datatype_parcellation(mesh);
headmodel = [];
if isfield(mesh,'pos')
  headmodel.pos = mesh.pos;
else
  ft_error('Vertex field is required!')
end

if isfield(mesh,'tet')
  headmodel.tet = mesh.tet;
elseif isfield(mesh,'hex')
  headmodel.hex = mesh.hex;
else
  ft_error('Connectivity information is required!')
end

if isfield(mesh,'tissue')
  headmodel.tissue = mesh.tissue;
else
  ft_error('No element indices declared!')
end

if isempty(conductivity)
  ft_error('No conductivity information!')
end

if length(conductivity) >= length(unique(headmodel.tissue))
  headmodel.cond = conductivity;
else
  ft_error('Wrong conductivity information!')
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


headmodel.stiff = sb_calc_stiff(headmodel);
headmodel.type  = 'simbio';
