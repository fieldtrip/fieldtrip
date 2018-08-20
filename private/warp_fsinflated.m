function [coord_inf] = warp_fsinflated(cfg, elec)

% WARP_FSINFLATED maps electrodes from FreeSurfer's pial surface to
% FreeSurfer's inflated brain.
%
% The configuration must contain the following options:
%   cfg.headshape      = string, filename containing subject headshape
%                      (e.g. <path to freesurfer/surf/lh.pial>)
%   cfg.fshome         = string, path to freesurfer
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2018, Richard Jimenez & Arjen Stolk
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
% $Id$

subj_pial = ft_read_headshape(cfg.headshape);
[PATHSTR, NAME] = fileparts(cfg.headshape); % lh or rh
if ~exist([PATHSTR filesep NAME '.inflated'],'file');
  ft_error([PATHSTR filesep NAME filesep 'inflated cannot be found'])
end
subj_inflated = ft_read_headshape([PATHSTR filesep NAME '.inflated']);

for e = 1:numel(elec.label)
  % subject space (3D surface): electrode pos -> pial vertex index
  dist = sqrt(sum(((subj_pial.pos - repmat(elec.elecpos(e,:), size(subj_pial.pos,1), 1)).^2),2));
  [~, minidx] = min(dist);
  
  % template space (3D surface): pial/inflated vertex index -> inflated brain pos
  coord_inf(e,:) = subj_inflated.pos(minidx,:);
  clear minidx
end
clear subj_pial subj_inflated
