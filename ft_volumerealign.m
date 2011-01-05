function [mri] = ft_volumerealign(cfg, mri)

% FT_VOLUMEREALIGN spatially aligns an anatomical MRI with head coordinates based on
% external fiducials or anatomical landmarks. This function does not change the
% volume itself, but adjusts the homogeneous transformation matrix that describes
% the coordinate system.
%
% This function only changes the coordinate system of an anatomical
% MRI, it does not change the MRI as such. For spatial normalisation
% (warping) of an MRI to a template brain you should use the
% FT_VOLUMENORMALISE function.
%
% Use as
%   [mri] = ft_volumerealign(cfg, mri)
% where mri is an anatomical volume (i.e. MRI) or a functional
% volume (i.e. source recunstruction that has been interpolated on
% an MRI).
%
% The configuration can contain the following options
%   cfg.clim           = [min max], scaling of the anatomy color (default
%                        is to adjust to the minimum and maximum)
%   cfg.method         = different methods for aligning the volume 
%                        'fiducial' realign the volume to the fiducials,
%                                     using 'ALS_CTF' convention, i.e.
%                                     the origin is exactly between lpa and rpa
%                                     the X-axis goes towards nas
%                                     the Y-axis goes approximately towards lpa, 
%                                       orthogonal to X and in the plane spanned
%                                       by the fiducials
%                                     the Z-axis goes approximately towards the vertex,
%                                       orthogonal to X and Y
%                        'landmark' realign the volume to anatomical landmarks,
%                                     using RAS_Tal convention, i.e.
%                                     the origin corresponds with the anterior commissure
%                                     the Y-axis is along the line from the posterior
%                                       commissure to the anterior commissure
%                                     the Z-axis is towards the vertex, in between the
%                                       hemispheres
%                                     the X-axis is orthogonal to the YZ-plane, 
%                                       positive to the right
%                        'interactive'     manually using graphical user interface
%
% For realigning to the fiducials, you should specify the position of the
% fiducials in voxel indices.
%   cfg.fiducial.nas  = [i j k], position of nasion
%   cfg.fiducial.lpa  = [i j k], position of LPA
%   cfg.fiducial.rpa  = [i j k], position of RPA
%
% For realigning to the landmarks, you should specify the position of the
% landmarks in voxel indices.
%   cfg.landmark.ac      = [i j k], position of anterior commissure
%   cfg.landmark.pc      = [i j k], position of posterior commissure
%   cfg.landmark.xzpoint = [i j k], point on the midsagittal-plane with positive Z-coordinate,
%                                     i.e. interhemispheric point above ac and pc
%
% See also FT_READ_MRI, FT_ELECTRODEREALIGN

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2006-2009, Robert Oostenveld
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

ft_defaults

cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial', 'fiducial'});
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'fiducial'),  cfg.fiducial = [];         end
if ~isfield(cfg, 'landmark'),  cfg.landmark = [];         end
if ~isfield(cfg, 'parameter'), cfg.parameter = 'anatomy'; end
if ~isfield(cfg, 'clim'),      cfg.clim      = [];        end
if ~isfield(cfg, 'inputfile'), cfg.inputfile = [];        end
if ~isfield(cfg, 'outputfile'),cfg.outputfile = [];       end

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    mri = loadvar(cfg.inputfile, 'mri');
  end
end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

if ~isfield(cfg, 'method')
  if ~isempty(cfg.fiducial)
    cfg.method = 'fiducial';
    basedonfid = 1;
    basedonmrk = 0;
  elseif ~isempty(cfg.landmark)
    cfg.method = 'landmark';
    basedonfid = 0;
    basedonmrk = 1;
  else
    cfg.method = 'interactive';
  end
end

if strcmp(cfg.method, 'interactive')
  basedonfid = 0;
  basedonmrk = 0;
end

% select the parameter that should be displayed
cfg.parameter = parameterselection(cfg.parameter, mri);
if iscell(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
end

switch cfg.method
  case 'fiducial'
    % do nothing
    
  case 'landmark'
    % do nothing

  case 'interactive'
    showcrosshair = true;
    dat = getsubfield(mri, cfg.parameter);
    nas = [];
    lpa = [];
    rpa = [];
    antcomm = [];
    pstcomm = [];
    xzpoint = [];
    x = 1:mri.dim(1);
    y = 1:mri.dim(2);
    z = 1:mri.dim(3);
    xc = round(mri.dim(1)/2);
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);
    
    while true % break when 'q' is pressed
      fprintf('============================================================\n');
      fprintf('click with mouse button to reslice the display to a new position\n');
      fprintf('press n/l/r on keyboard to record the current position as fiducial location\n');
      fprintf('press a/p/z on keyboard to record the current position as anatomical landmark\n');
      fprintf('press i,j,k or I,J,K on the keyboard to increment or decrement the slice number by one\n');
      fprintf('press c or C on the keyboard to show or hide the crosshair\n');
      fprintf('press q on keyboard to quit interactive mode\n');
      xc = round(xc); xc = max(1,xc); xc = min(mri.dim(1),xc);
      yc = round(yc); yc = max(1,yc); yc = min(mri.dim(2),yc);
      zc = round(zc); zc = max(1,zc); zc = min(mri.dim(3),zc);
      volplot(x, y, z, dat, [xc yc zc], cfg.clim, showcrosshair);
      drawnow;
      try, [d1, d2, key] = ginput(1); catch, key='q'; end
      if key=='q'
        break;
      elseif key=='l'
        lpa = [xc yc zc];
      elseif key=='r'
        rpa = [xc yc zc];
      elseif key=='n'
        nas = [xc yc zc];
      elseif key=='a'
        antcomm = [xc yc zc];
      elseif key=='p'
        pstcomm = [xc yc zc];
      elseif key=='z'
        xzpoint = [xc yc zc];
      elseif key=='c'
        showcrosshair = true;
      elseif key=='C'
        showcrosshair = false;
      elseif key=='i'
        xc = xc+1;
      elseif key=='I'
        xc = xc-1;
      elseif key=='j'
        yc = yc+1;
      elseif key=='J'
        yc = yc-1;
      elseif key=='k'
        zc = zc+1;
      elseif key=='K'
        zc = zc-1;
      else
        % update the view to a new position
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        switch l1,
          case 'i'
            xc = d1;
          case 'j'
            yc = d1;
          case 'k'
            zc = d1;
        end
        switch l2,
          case 'i'
            xc = d2;
          case 'j'
            yc = d2;
          case 'k'
            zc = d2;
        end
      end
      
      fprintf('cur_voxel = [%f %f %f], cur_head = [%f %f %f]\n', [xc yc zc], warp_apply(mri.transform, [xc yc zc]));
      if ~isempty(nas),
        fprintf('nas_voxel = [%f %f %f], nas_head = [%f %f %f]\n', nas, warp_apply(mri.transform, nas));
      else
        fprintf('nas = undefined\n');
      end
      if ~isempty(lpa),
        fprintf('lpa_voxel = [%f %f %f], lpa_head = [%f %f %f]\n', lpa, warp_apply(mri.transform, lpa));
      else
        fprintf('lpa = undefined\n');
      end
      if ~isempty(rpa),
        fprintf('rpa_voxel = [%f %f %f], rpa_head = [%f %f %f]\n', rpa, warp_apply(mri.transform, rpa));
      else
        fprintf('rpa = undefined\n');
      end
    end % while true
    
    cfg.fiducial.nas = nas;
    cfg.fiducial.lpa = lpa;
    cfg.fiducial.rpa = rpa;
    
    cfg.landmark.ac     = antcomm;
    cfg.landmark.pc     = pstcomm;
    cfg.landmark.xzpoint = xzpoint;

    if ~isempty(nas) && ~isempty(lpa) && ~isempty(rpa)
      basedonfid = 1;
    end

    if ~isempty(antcomm) && ~isempty(pstcomm) && ~isempty(xzpoint)
      basedonmrk = 1;
    end
  otherwise
    error('unsupported method');
end

if basedonfid && basedonmrk
  basedonmrk = 0;
  warning('both fiducials and anatomical landmarks have been defined interactively: using the fiducials for realignment');
end

if basedonfid
  % the fiducial locations are now specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  nas_head = warp_apply(mri.transform, cfg.fiducial.nas);
  lpa_head = warp_apply(mri.transform, cfg.fiducial.lpa);
  rpa_head = warp_apply(mri.transform, cfg.fiducial.rpa);

  % compute the homogenous transformation matrix describing the new coordinate system
  [realign, coordsys] = headcoordinates(nas_head, lpa_head, rpa_head);

elseif basedonmrk
  % the fiducial locations are now specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  ac     = warp_apply(mri.transform, cfg.landmark.ac);
  pc     = warp_apply(mri.transform, cfg.landmark.pc);
  xzpoint= warp_apply(mri.transform, cfg.landmark.xzpoint);

  % compute the homogenous transformation matrix describing the new coordinate system
  [realign, coordsys] = headcoordinates(ac, pc, xzpoint, 'RAS_TAL');

end

% combine the additional transformation with the original one
mri.transform = realign * mri.transform;
mri.coordsys  = coordsys;

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();

% remember the configuration
mri.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'mri', mri); % use the variable name "data" in the output file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to show three orthogonal slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function volplot(x, y, z, dat, c, cscale, showcrosshair)
xi = c(1);
yi = c(2);
zi = c(3);

% manual color scaling of anatomy data is usefull in case of some pixel noise
if nargin<6 || isempty(cscale)
  cmin = min(dat(:));
  cmax = max(dat(:));
else
  cmin = cscale(1);
  cmax = cscale(2);
end

h1 = subplot(2,2,1);
h2 = subplot(2,2,2);
h3 = subplot(2,2,3);

subplot(h1);
imagesc(x, z, squeeze(dat(:,yi,:))'); set(gca, 'ydir', 'normal')
axis equal; axis tight;
xlabel('i'); ylabel('k');
caxis([cmin cmax]);
if showcrosshair
  crosshair([x(xi) z(zi)], 'color', 'yellow');
end

subplot(h2);
imagesc(y, z, squeeze(dat(xi,:,:))'); set(gca, 'ydir', 'normal')
axis equal; axis tight;
xlabel('j'); ylabel('k');
caxis([cmin cmax]);
if showcrosshair
  crosshair([y(yi) z(zi)], 'color', 'yellow');
end

subplot(h3);
imagesc(x, y, squeeze(dat(:,:,zi))'); set(gca, 'ydir', 'normal')
axis equal; axis tight;
xlabel('i'); ylabel('j');
caxis([cmin cmax]);
if showcrosshair
  crosshair([x(xi) y(yi)], 'color', 'yellow');
end

colormap gray
