function [coord_norm] = warp_fsaverage_sym(cfg, elec)

% WARP_FSAVERAGE_SYM maps left or right hemisphere electrodes onto 
% FreeSurfer's fsaverage_sym's left hemisphere. To perform this mapping, 
% you first need to have processed the subject's MRI with FreeSurfer's 
% recon-all functionality and additionaly have registered the subject's resulting 
% surfaces to freesurfer fsaverage_sym template using surfreg as described 
% in section 1.2 of https://surfer.nmr.mgh.harvard.edu/fswiki/Xhemi
%
% The configuration must contain the following options
%   cfg.headshape      = string, filename containing subject headshape
%                      (e.g. <path to freesurfer/surf/lh.pial>)
%   cfg.fshome         = string, path to freesurfer
%
% See also FT_ELECTRODEREALIGN, WARP_FSAVERAGE

% Copyright (C) 2019, Arjen Stolk
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
if strcmp(NAME, 'lh')
  subj_reg = ft_read_headshape([PATHSTR filesep 'lh.fsaverage_sym.sphere.reg']);
elseif strcmp(NAME, 'rh')
  subj_reg = ft_read_headshape([PATHSTR(1:strfind(PATHSTR, [filesep 'surf'])-1) filesep 'xhemi' filesep 'surf' filesep 'lh.fsaverage_sym.sphere.reg']);
end
if ~isfolder([cfg.fshome filesep 'subjects' filesep 'fsaverage_sym']) || ~isfolder([PATHSTR(1:strfind(PATHSTR, [filesep 'surf'])-1) filesep 'xhemi'])
  ft_error(['fsaverage_sym and/or xhemi folders cannot be found'])
end
fsavg_pial = ft_read_headshape([cfg.fshome filesep 'subjects' filesep 'fsaverage_sym' filesep 'surf' filesep 'lh.pial']);
fsavg_reg = ft_read_headshape([cfg.fshome filesep 'subjects' filesep 'fsaverage_sym' filesep 'surf' filesep 'lh.sphere.reg']); % always map onto the left hemi

for e = 1:numel(elec.label)
  % subject space (3D surface): electrode pos -> vertex index
  dist = sqrt(sum(((subj_pial.pos - repmat(elec.elecpos(e,:), size(subj_pial.pos,1), 1)).^2),2));
  [dum, minidx] = min(dist);
  
  % intersubject space (3D sphere): vertex index -> vertex pos -> template vertex index
  dist2 = sqrt(sum(((fsavg_reg.pos - repmat(subj_reg.pos(minidx,:), size(fsavg_reg.pos,1), 1)).^2),2));
  [dum, minidx2] = min(dist2);
  clear minidx
  
  % template space (3D surface): template vertex index -> template electrode pos
  coord_norm(e,:) = fsavg_pial.pos(minidx2,:);
  clear minidx2
end
clear subj_pial subj_reg fsavg_pial fsavg_reg
