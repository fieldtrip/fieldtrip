function vol = ft_headmodel_halfspace(location, orientation, varargin)

% FT_HEADMODEL_HALFSPACE creates an EEG volume conduction model that
% is described with an infinite conductive halfspace. You can think
% of this as a plane with on one side a infinite mass of conductive
% material (e.g. water) and on the other side non-conductive material
% (e.g. air).
% 
% Use as
%   vol = ft_headmodel_halfspace(location, orientation, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   location    = 1x3 vector specifying a point on the plane 
%   orientation = 1x3 vector specifying the direction orthogonal to the plane
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% the description of this volume conduction model only consists of the desctiption of the plane.
vol = [];
vol.pos = location;    % a point that lies in the plane that separates the conductive tissue from the air
vol.ori = orientation; % a unit vector pointing towards the air
vol.ori = vol.ori/norm(vol.ori);
vol.type = 'halfspace';

