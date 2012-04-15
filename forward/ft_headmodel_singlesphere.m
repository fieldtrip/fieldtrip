function vol = ft_headmodel_singlesphere(geometry, varargin)

% FT_HEADMODEL_SINGLESPHERE creates a volume conduction model of the
% head by fitting a spherical model to a set of points that describe
% the head surface.
%
% For MEG this implements Cuffin BN, Cohen D.  "Magnetic fields of
% a dipole in special volume conductor shapes" IEEE Trans Biomed Eng.
% 1977 Jul;24(4):372-81.
%
% Use as
%   vol = ft_headmodel_singlesphere(pnt, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   conductivity     = number, conductivity of the sphere
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% FIXME document the EEG case

% get the optional arguments
conductivity = ft_getopt(varargin, 'conductivity',1);
unit         = ft_getopt(varargin,'unit');

if length(conductivity)~=1
  error('the conductivity should be a single number')
end

% start with an empty volume conductor
vol = [];

if ~isempty(unit)
  % use the user-specified units for the output
  vol.unit = geometry.unit;
elseif isfield(geometry, 'unit')
  % copy the geometrical units into he volume conductor
  vol.unit = geometry.unit;
end

if isnumeric(geometry) && size(geometry,2)==3
  % assume that it is a Nx3 array with vertices
elseif isstruct(geometry) && isfield(geometry,'pnt') && numel(geometry)==1
  % get the points from the triangulated surface
  geometry = geometry.pnt;
else
  error('the input geometry should be a set of points or a single triangulated surface')
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(geometry);

vol.r = single_r;
vol.o = single_o;
vol.c = conductivity;
vol.type = 'singlesphere';
vol      = ft_convert_units(vol); % ensure the object to have a unit

