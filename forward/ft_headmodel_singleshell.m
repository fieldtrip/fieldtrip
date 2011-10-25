function vol = ft_headmodel_singleshell(geometry, varargin)

% FT_HEADMODEL_SINGLESHELL creates a volume conduction model of the
% head for MEG based on a realistic shaped surface of the inside of
% the skull.
% 
% The method implemented in this function allows for a simple and
% fast method for the MEG forward calculation for one shell of arbitrary
% shape, based on a correction of the lead field for a spherical
% volume conductor by a superposition of basis functions, gradients
% of harmonic functions constructed from spherical harmonics.
% 
% This function implements
%   G. Nolte, "The magnetic lead field theorem in the quasi-static
%   approximation and its use for magnetoencephalography forward calculation
%   in realistic volume conductors", Phys Med Biol. 2003 Nov 21;48(22):3637-52.
% 
% Use as
%   vol = ft_headmodel_singleshell(geom, ...)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% if it contains more than 1 shell it retunrs an error
if isfield(geometry,'pnt')
  if numel(geometry)>1
    error('no more than 1 shell at a time is allowed')
  end
else
  error('the input should be a boundary')
end

vol      = [];
vol.bnd  = geometry;
vol.type = 'nolte';
if ~isfield(vol, 'unit')
  vol = ft_convert_units(vol);
end
