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

% this is an easy one
vol = [];
vol.type = 'infinite';

