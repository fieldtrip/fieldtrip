function [coord_norm] = warp_fsaverage(cfg, elec)

% WARP_FSAVERAGE maps electrodes onto FreeSurfer's fsaverage brain.
% This surface-based registration technique solely considers the curvature
% patterns of the cortex and thus can be used for the spatial normalization
% of electrodes located on or near the cortical surface. To perform
% surface-based normalization, you first need to process the subject's MRI
% with FreeSurfer's recon-all functionality.
%
% The configuration must contain the following options
%   cfg.headshape      = string, filename containing subject headshape 
%                      (e.g. <path to freesurfer/surf/lh.pial>)
%   cfg.fshome         = string, path to freesurfer
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2017, Arjen Stolk
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
subj_reg = ft_read_headshape([PATHSTR filesep NAME '.sphere.reg']);
if ~isfolder([cfg.fshome filesep 'subjects' filesep 'fsaverage' filesep 'surf'])
  ft_error(['freesurfer dir ' cfg.fshome filesep 'subjects' filesep 'fsaverage' filesep 'surf cannot be found'])
end
fsavg_pial = ft_read_headshape([cfg.fshome filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep NAME '.pial']);
fsavg_reg = ft_read_headshape([cfg.fshome filesep 'subjects' filesep 'fsaverage' filesep 'surf' filesep NAME '.sphere.reg']);

for e = 1:numel(elec.label)
  % subject space (3D surface): electrode pos -> vertex index
  dist = sqrt(sum(((subj_pial.pos - repmat(elec.elecpos(e,:), size(subj_pial.pos,1), 1)).^2),2));
  [~, minidx] = min(dist);
  
  % intersubject space (2D sphere): vertex index -> vertex pos -> template vertex index
  dist2 = sqrt(sum(((fsavg_reg.pos - repmat(subj_reg.pos(minidx,:), size(fsavg_reg.pos,1), 1)).^2),2));
  [~, minidx2] = min(dist2);
  clear minidx
  
  % template space (3D surface): template vertex index -> template electrode pos
  coord_norm(e,:) = fsavg_pial.pos(minidx2,:);
  clear minidx2
end
clear subj_pial subj_reg fsavg_pial fsavg_reg
