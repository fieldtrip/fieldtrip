function vol = ft_headmodel_simbio(geom, varargin)

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
%   vol = ft_headmodel_simbio(geom,'conductivity', conductivities, ...)
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
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% $Id$

ft_hastoolbox('simbio', 1);

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');

% start with an empty volume conductor
geom = ft_datatype_parcellation(geom);
vol = [];
if isfield(geom,'pos')
  vol.pos = geom.pos;
else
  error('Vertex field is required!')
end

if isfield(geom,'tet')
  vol.tet = geom.tet;
elseif isfield(geom,'hex')
  vol.hex = geom.hex;
else
  error('Connectivity information is required!')
end

if isfield(geom,'tissue')
  vol.tissue = geom.tissue;
else
  error('No element indices declared!')
end

if isempty(conductivity)
  error('No conductivity information!')
end

if length(conductivity) >= length(unique(vol.tissue))
  vol.cond = conductivity;
else
  error('Wrong conductivity information!')
end

if ~isfield(geom,'tissuelabel')
        numlabels = size(unique(geom.tissue),1);
    vol.tissuelabel = {};
            ulabel = unique(labels);
    for i = 1:numlabels
        shape.tissuelabel{i} = num2str(ulabel(i));
    end
else
    vol.tissuelabel = geom.tissuelabel;
end
    
vol.stiff = sb_calc_stiff(vol);
vol.type = 'simbio';

end


