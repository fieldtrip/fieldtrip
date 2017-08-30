function [surface_original] = ft_surfacerealign(cfg, surface_original)

% FT_SURFACEREALIGN realigns surface
%
% FIDUCIAL - You can apply a rigid body realignment based on three fiducial
% locations. After realigning, the fiducials in the input surface
% (typically nose, left and right ear) are along the same axes as the
% fiducials in the template surface set.
%
%   cfg.method         = string representing the method for aligning the surface
%                        'fiducial'        realign using three fiducials
%                        (e.g. NAS, LPA and RPA)
%
% If you want to realign the surface using fiducials, the target and the
% objective have to contain the three fiducials which relate , e.g.
%   cfg.target.elecpos(1,:)    = [110 0 0]     % location of the nose
%   cfg.target.elecpos(2,:)    = [0  90 0]     % location of the left ear
%   cfg.target.elecpos(3,:)    = [0 -90 0]     % location of the right ear
%   cfg.target.label       = {'NAS', 'LPA', 'RPA'}
%   cfg.objective.elecpos(1,:) = [0 -110 0]      % location of the nose
%   cfg.objective.elecpos(2,:) = [90   0 0]      % location of the left ear
%   cfg.objective.elecpos(3,:) = [-90  0 0]      % location of the right ear
%   cfg.objective.label    = {'NAS', 'LPA', 'RPA'}
%
% See also FT_MESHREALIGN

% Copyright (C) 2016, Simon Homoelle
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

% DEPRECATED by roboos on 22 August 2017
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1830
% support for this functionality can be removed mid 2018
warning('FT_SURFACEREALIGN is deprecated, please use FT_MESHREALIGN instead.')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision              = '$Id$';
ft_nargin                = nargin;
ft_nargout               = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    surface_original
ft_preamble provenance surface_original
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

usefiducial              = isfield(cfg, 'target') & isfield(cfg, 'objective') & strcmp(cfg.method,'fiducial');

if usefiducial
  cfg.elec                 = cfg.objective;
  %use ft_electroderealign for realign the surface
  surface_realigned           = ft_electroderealign(cfg);
  %use transformation obtained by ft_electroderealign
  transform                = [surface_original.pos, ones(length(surface_original.pos),1)]*surface_realigned.homogeneous';
  surface_original.pos        = transform(:,1:3);
  surface_original.cfg        = surface_realigned.cfg;
else
  ft_error('Cannot perform ft_surfacerealign. Please read help ft_surfacerealign, and check your cfg')
end

