function headmodel = ft_headmodel_simbio(geom, varargin)

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
%   headmodel = ft_headmodel_simbio(geom,'conductivity', conductivities, ...)
%
% The geom is given as a volumetric mesh, using ft_datatype_parcellation
%   geom.pos = vertex positions
%   geom.tet/geom.hex = list of volume elements
%   geom.tissue = tissue assignment for elements
%   geom.tissuelabel = labels correspondig to tissues
%
% Required input arguments should be specified in key-value pairs and have
% to include
%   conductivity   = vector containing tissue conductivities using ordered
%                    corresponding to geom.tissuelabel
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
geom = ft_datatype_parcellation(geom);
headmodel = [];
if isfield(geom,'pos')
  headmodel.pos = geom.pos;
else
  error('Vertex field is required!')
end

if isfield(geom,'tet')
  headmodel.tet = geom.tet;
elseif isfield(geom,'hex')
  headmodel.hex = geom.hex;
else
  error('Connectivity information is required!')
end

if isfield(geom,'tissue')
  headmodel.tissue = geom.tissue;
else
  error('No element indices declared!')
end

if isempty(conductivity)
  error('No conductivity information!')
end

if length(conductivity) >= length(unique(headmodel.tissue))
  headmodel.cond = conductivity;
else
  error('Wrong conductivity information!')
end

if ~isfield(geom,'tissuelabel')
  numlabels = size(unique(geom.tissue),1);
  headmodel.tissuelabel = {};
  ulabel = unique(labels);
  for i = 1:numlabels
    headmodel.tissuelabel{i} = num2str(ulabel(i));
  end
else
  headmodel.tissuelabel = geom.tissuelabel;
end

headmodel.stiff = sb_calc_stiff(headmodel);
headmodel.type = 'simbio';

end
