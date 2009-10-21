function [vol, cfg] = prepare_localspheres(cfg, mri)

% PREPARE_LOCALSPHERES creates a MEG volume conductor model with a sphere
% for every sensor. You can also use it to create a single sphere
% model that is fitted to the MRI or to the head shape points.
%
% Use as
%   [vol, cfg] = prepare_localspheres(cfg, seg), or
%   [vol, cfg] = prepare_localspheres(cfg, mri), or
%   [vol, cfg] = prepare_localspheres(cfg)
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
%   cfg.mriunits    = 'mm' or 'cm' (default = 'mm')
%   cfg.sourceunits = 'mm' or 'cm' (default = 'cm')
%   cfg.threshold   = 0.5, relative to the maximum value in the segmentation
%
% This function implements
%   Huang MX, Mosher JC, Leahy RM.
%   A sensor-weighted overlapping-sphere head model and exhaustive head model comparison for MEG
%   Phys Med Biol. 1999 Feb;44(2):423-40

% TODO cfg.spheremesh  should be renamed consistently with other mesh generation cfgs
% TODO shape should contain pnt as subfield and not be equal to pnt (for consistency with other use of shape)
%
% Undocumented local options:
% cfg.spheremesh, number of points that is placed on the brain surface (default 4000)
% cfg.maxradius

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: prepare_localspheres.m,v $
% Revision 1.29  2009/07/29 06:40:41  roboos
% updated plotting functions
%
% Revision 1.28  2009/07/16 09:17:17  crimic
% link to prepare_mesh.m function
%
% Revision 1.27  2009/05/29 10:47:19  roboos
% only convert to struct in case headshape is specified
%
% Revision 1.26  2009/05/25 08:04:40  roboos
% fixed the name of the "headshape" variable
% ensure that cfg.headshape is a structure and not a config object
%
% Revision 1.25  2009/05/18 13:59:14  vlalit
% typo fix
%
% Revision 1.24  2009/05/14 19:21:36  roboos
% consistent handling of cfg.headshape in code and documentation
%
% Revision 1.23  2009/04/01 12:29:26  roboos
% added checkconfig
%
% Revision 1.22  2009/04/01 12:17:24  roboos
% use Taubin's method instead of optimization toolbox for fitting the sphere (thanks to Jean and Guillaume)
%
% Revision 1.21  2009/01/07 14:28:51  roboos
% do not open new figure but start with "clf"
% removed typicalx from optimization
%
% Revision 1.20  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.19  2008/08/13 21:02:20  roboos
% use general read_headshape instead of specific subfunctions
%
% Revision 1.18  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.17  2007/08/06 09:20:14  roboos
% added support for bti_hs
%
% Revision 1.16  2007/07/26 08:02:10  roboos
% remove double vertices in headshape
% changed some indentation
%
% Revision 1.15  2006/10/09 15:25:00  roboos
% added option to fit only a single sphere to the MRI or headshape points
%
% Revision 1.14  2006/08/16 10:52:12  marsie
% fixed use of spm_smooth()
%
% Revision 1.13  2006/08/01 10:32:01  marsie
% fixed bug in using spm_smooth
%
% Revision 1.12  2006/07/27 08:29:38  roboos
% use spm_smooth instead of spm_conv, updated documentation
%
% Revision 1.11  2006/06/08 07:50:52  roboos
% updated the conversion between the source and MRI units (support mm,cm,dm,m for both)
%
% Revision 1.10  2006/06/07 09:34:19  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.9  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.8  2006/04/18 19:04:35  roboos
% changed hard-coded 10000 vertex points for brain surface into cfg.spheremesh with default value of 4000
%
% Revision 1.7  2006/04/05 15:07:37  roboos
% return the configuration as second argument, store the headshape points in the cfg
%
% Revision 1.6  2006/04/03 10:47:13  roboos
% pre-allocate memory to hold the result, ensure that the dimensions are correct (thanks to Hanneke)
%
% Revision 1.5  2006/03/21 09:45:07  roboos
% some irrelevant changes in the code flow to deal with presence/absence of segmented MRI input
%
% Revision 1.4  2006/03/09 08:18:34  roboos
% added code to use a segmented anatomical MRI, i.e. smooth and do edge detection
% with options similar to prepare_dipole_grid in gray matter
% updated help
%
% Revision 1.3  2006/02/06 21:22:28  roboos
% removed obvious dependency on CTF for reading hdr.grad
%
% Revision 1.2  2006/01/25 10:08:27  roboos
% added the complete implementation of this function, sofar it had been empty
% the implementation is based on code from jansch, which again was partially based on roboos' code
%
% Revision 1.1  2005/11/03 11:15:45  roboos
% new implementation
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'radius'),        cfg.radius = 8.5;        end
if ~isfield(cfg, 'maxradius'),     cfg.maxradius = 20;      end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 5;        end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';    end
if ~isfield(cfg, 'smooth');        cfg.smooth    = 5;       end % in voxels
if ~isfield(cfg, 'mriunits');      cfg.mriunits = 'mm';     end
if ~isfield(cfg, 'sourceunits'),   cfg.sourceunits = 'cm';  end
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 4000;   end
if ~isfield(cfg, 'singlesphere'),  cfg.singlesphere = 'no'; end
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];      end

% construct the geometry of the headshape using a single boundary
if nargin==1
  headshape = prepare_mesh(cfg);
else
  headshape = prepare_mesh(cfg, mri);
end

% read the gradiometer definition from file or copy it from the configuration
if isfield(cfg, 'gradfile')
  grad = read_sens(cfg.gradfile);
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
  plot_sens(grad);
  plot_mesh(headshape, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
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
  allpnt  = grad.pnt(coilsel, :);	% position of all coils belonging to this channel
  allori  = grad.ori(coilsel, :);	% orientation of all coils belonging to this channel
  
  if strcmp(cfg.feedback, 'yes')
    cla
    plot3(grad.pnt(:,1), grad.pnt(:,2), grad.pnt(:,3), 'b.');	% all coils
    plot3(allpnt(:,1), allpnt(:,2), allpnt(:,3), 'r*');		% this channel in red
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
    plot_mesh(headshape.pnt(shapesel,:), 'vertexcolor', 'g');
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

vol.type = 'multisphere';

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

