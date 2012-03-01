function vol = ft_headmodel_infinite(varargin)

% FT_HEADMODEL_INFINITE returns an infinitely large homogenous
% volume conduction model. For EEG the volume conductor can be used
% to compute the leadfield of electric current dipoles, for MEG it
% can be used for computing the leadfield of magnmetic dipoles.
% 
% Use as
%   vol = ft_headmodel_infinite;
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

model = ft_getopt(varargin, 'sourcemodel', 'dipole');

% this is an easy one
vol = [];

if strcmp(model,'monopole')
  vol.type = 'infinite_monopole';  
else
  vol.type = 'infinite';
end
