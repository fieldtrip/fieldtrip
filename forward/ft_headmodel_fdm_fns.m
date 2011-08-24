function vol = ft_headmodel_fdm_fns(varargin)

% FT_HEADMODEL_FDM_FNS creates the volume conduction structure to be used 
% in the FNS forward solver.
%
% Use as
%   vol = ft_headmodel_fdm_fns(varargin)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = matrix C [9XN tissue types]
%
% The conductivity matrix C is stored in a 9xN table, where N is the number of tissues and a 
% 3x3 tensor conductivity matrix is stored in each row in column oriented format
% 
% Standard default values for conductivity matrix C are derived from 
% Saleheen HI, Ng KT. New finite difference formulations for general
% inhomogeneous anisotropic bioelectric problems. IEEE Trans Biomed Eng.
% 1997
% 
% Additional documentation available at:
% http://hunghienvn.nmsu.edu/wiki/index.php/FNS
% 
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% get the optional arguments
conductivity = keyval('conductivity', varargin);

if isempty(conductivity)
  conductivity = fns_contable_write;
end

% start with an empty volume conductor
vol = [];
vol.cond = conductivity;
vol.type = 'fns';
