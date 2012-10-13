function [vol, cfg] = ft_prepare_localspheres(cfg, mri)

% FT_PREPARE_LOCALSPHERES creates a MEG volume conductor model with a sphere
% for every sensor. You can also use it to create a single sphere
% model that is fitted to the MRI or to the head shape points.
%
% Use as
%   vol = ft_prepare_localspheres(cfg, seg), or
%   vol = ft_prepare_localspheres(cfg, mri), or
%   vol = ft_prepare_localspheres(cfg)
%
% The input configuration should contain
%   cfg.grad         = structure with gradiometer definition, or
%   cfg.gradfile     = filename containing gradiometer definition
%   cfg.radius       = number, which points to select for each channel (default = 7 cm)
%   cfg.baseline     = number, baseline of axial/planar gradiometer (default = 5 cm)
%   cfg.feedback     = 'yes' or 'no' (default = 'yes')
%   cfg.singlesphere = 'yes' or 'no', fit only a single sphere (default = 'no')
%   cfg.headshape    = a filename containing headshape, a structure containing a
%                      single triangulated boundary, or a Nx3 matrix with surface
%                      points
%
% The following options are relevant if you use a segmented MRI
%   cfg.smooth      = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%   cfg.sourceunits = 'mm' or 'cm' (default = 'cm')
%   cfg.threshold   = 0.5, relative to the maximum value in the segmentation
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% This function implements
%   Huang MX, Mosher JC, Leahy RM.
%   A sensor-weighted overlapping-sphere head model and exhaustive head model comparison for MEG
%   Phys Med Biol. 1999 Feb;44(2):423-40

% TODO cfg.spheremesh  should be renamed consistently with other mesh generation cfgs
% TODO shape should contain pnt as subfield and not be equal to pnt (for consistency with other use of shape)
%
% Undocumented local options:
% cfg.spheremesh = number of points that is placed on the brain surface (default 4000)
% cfg.maxradius
%
% See also FT_PREPARE_CONCENTRICSPHERES, FT_PREPARE_BEMMODEL,
% FT_PREPARE_SINGLESHELL, FT_PREPARE_LEADFIELD, FT_PREPARE_MESH,
% FT_PREPARE_MESH_NEW

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar mri

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'spheremesh', 'numvertices'}); 

% set the defaults
if ~isfield(cfg, 'radius'),        cfg.radius = 8.5;        end
if ~isfield(cfg, 'maxradius'),     cfg.maxradius = 20;      end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 5;        end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';    end
if ~isfield(cfg, 'smooth');        cfg.smooth    = 5;       end % in voxels
if ~isfield(cfg, 'sourceunits'),   cfg.sourceunits = 'cm';  end
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'numvertices'),   cfg.numvertices = [];    end
if ~isfield(cfg, 'singlesphere'),  cfg.singlesphere = 'no'; end
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];      end

hasmri = exist('mri', 'var'); % note that nargin will not work in case of cfg.inputfile
if hasmri
  headshape = ft_prepare_mesh(cfg, mri);
elseif isfield(cfg,'headshape') && nargin == 1 
  if ischar(cfg.headshape)
    headshape = ft_read_headshape(cfg.headshape);
  else
    headshape = cfg.headshape;
  end
else
  error('no head shape available')
end

% read the gradiometer definition from file or copy it from the configuration
if isfield(cfg, 'gradfile')
  grad = ft_read_sens(cfg.gradfile);
else
  grad = cfg.grad;
end

Nshape = size(headshape.pnt,1);
Nchan  = size(grad.tra, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
end

% plot all channels and headshape points
if strcmp(cfg.feedback, 'yes')
  cla
  ft_plot_sens(grad);
  ft_plot_mesh(headshape, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
  drawnow
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(headshape.pnt);
fprintf('single sphere,   %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', Nshape, single_o(1), single_o(2), single_o(3), single_r);

vol = [];

if strcmp(cfg.singlesphere, 'yes')
  % only return a single sphere
  vol.r = single_r;
  vol.o = single_o;
  return;
end

% start with an empty structure that will hold the results
vol.r = zeros(Nchan,1);    % radius of every sphere
vol.o = zeros(Nchan,3);    % origin of every sphere
vol.label = cell(Nchan,1); % corresponding gradiometer channel label for every sphere

for chan=1:Nchan
  coilsel = find(grad.tra(chan,:)~=0);
  allpnt  = grad.coilpos(coilsel, :);   % position of all coils belonging to this channel
  allori  = grad.coilori(coilsel, :);   % orientation of all coils belonging to this channel
  
  if strcmp(cfg.feedback, 'yes')
    cla
    plot3(grad.coilpos(:,1), grad.coilpos(:,2), grad.coilpos(:,3), 'b.');   % all coils
    plot3(allpnt(:,1), allpnt(:,2), allpnt(:,3), 'r*');     % this channel in red
  end
  
  % determine the average position and orientation of this channel
  thispnt = mean(allpnt,1);
  [u, s, v] = svd(allori);
  thisori = v(:,1)';
  if dot(thispnt,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end
  
  % compute the distance from every coil along this channels orientation
  dist = zeros(size(coilsel));
  for i=1:length(coilsel)
    dist(i) = dot((allpnt(i,:)-thispnt), thisori);
  end
  
  [m, i] = min(dist);
  % check whether the minimum difference is larger than a typical distance
  if abs(m)>(cfg.baseline/4)
    % replace the position of this channel by the coil that is the closest to the head (axial gradiometer)
    % except when the center of the channel is approximately just as good (planar gradiometer)
    thispnt = allpnt(i,:);
  end
  
  % find the headshape points that are close to this channel
  dist = sqrt(sum((headshape.pnt-repmat(thispnt,Nshape,1)).^2, 2));
  shapesel = find(dist<cfg.radius);
  if strcmp(cfg.feedback, 'yes')
    ft_plot_mesh(headshape.pnt(shapesel,:), 'vertexcolor', 'g');
    drawnow
  end
  
  % fit a sphere to these headshape points
  if length(shapesel)>10
    [o, r] = fitsphere(headshape.pnt(shapesel,:));
    fprintf('channel = %s, %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', grad.label{chan}, length(shapesel), o(1), o(2), o(3), r);
  else
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end
  
  if r > cfg.maxradius
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end
  
  % add this sphere to the volume conductor
  vol.o(chan,:)   = o;
  vol.r(chan)     = r;
  vol.label{chan} = grad.label{chan};
end % for all channels

vol.type = 'localspheres';

% ensure that the geometrical units are specified
vol = ft_convert_units(vol);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance
if hasmri
  ft_postamble previous mri
end
ft_postamble history vol

