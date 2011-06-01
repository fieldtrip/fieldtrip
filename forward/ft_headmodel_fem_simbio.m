function vol = ft_headmodel_fem_simbio(varargin)
% FT_HEADMODEL_FEM_SIMBIO reads a volume conduction model from a Vista .v
% file
%
% Vista is a software ...
% 
% Use as
%   vol = ft_headmodel_fem_simbio(filename)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% read the headmodel from file
% this works for Vista version x.x
vol = ft_read_vol(filename);