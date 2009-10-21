function [mri] = volumerealign(cfg, mri);

% VOLUMEREALIGN spatially aligns an anatomical MRI with head coordinates based on
% external fiducials. This function does not change the volume
% itself, but adjusts the homogenous transformation matrix that
% describes the coordinate system.
%
% This function only changes the coordinate system of an anatomical
% MRI, it does not change the MRI as such. For spatial normalisation
% (warping) of an MRI to a template brain you should use the
% VOLUMENORMALISE function.
%
% Use as
%   [mri] = volumerealign(cfg, mri)
% where mri is an anatomical volume (i.e. MRI) or a functional
% volume (i.e. source recunstruction that has been interpolated on
% an MRI).
%
% The configuration can contain the following options
%   cfg.clim           = [min max], scaling of the anatomy color (default
%                        is to adjust to the minimum and maximum)
%   cfg.method         = different methods for aligning the electrodes
%                        'realignfiducial' realign the volume to the fiducials
%                        'interactive'     manually using graphical user interface
%
% For realigning to the fiducials, you should specify the position of the
% fiducials in voxel indices.
%   cfg.fiducial.nas  = [i j k], position of nasion
%   cfg.fiducial.lpa  = [i j k], position of LPA
%   cfg.fiducial.rpa  = [i j k], position of RPA
%
% By specifying the fiducial coordinates either in the cfg or interactively, the
% anatomical MRI volume is realigned according to the folowing convention:
% - the origin is exactly between LPA and RPA
% - the X-axis goes towards NAS
% - the Y-axis goes approximately towards LPA, orthogonal to X and in the plane spanned by the fiducials
% - the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
%
% See also READ_MRI, ELECTRODEREALIGN

% Copyright (C) 2006-2009, Robert Oostenveld
%
% $Log: volumerealign.m,v $
% Revision 1.12  2009/07/31 13:43:36  jansch
% now really fixed a bug (unlike last time)
%
% Revision 1.11  2009/07/30 14:22:00  jansch
% fixed bug in input arguments for volplot
%
% Revision 1.10  2009/07/29 13:53:49  roboos
% added cfg.clim for color scaling, thanks to Hanneke
%
% Revision 2    2009/07/29 13:00:00 hanvdij
% Added colorscaling option in the volplot function at the end of this
% script.
%
% Revision 1.9  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.8  2008/11/21 13:56:12  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.7  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.6  2007/05/02 15:22:37  roboos
% cfg.parameter should never be a cell
%
% Revision 1.5  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.4  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.3  2006/10/10 16:21:10  roboos
% fixed bug in default setting for method
%
% Revision 1.2  2006/10/10 13:38:31  roboos
% added some help, thanks to Till
%
% Revision 1.1  2006/10/10 10:25:59  roboos
% new impementation
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
mri = checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'fiducial'),  cfg.fiducial = [];         end
if ~isfield(cfg, 'parameter'), cfg.parameter = 'anatomy'; end
if ~isfield(cfg, 'clim'),      cfg.clim      = [];        end

if ~isfield(cfg, 'method')
  if ~isempty(cfg.fiducial)
    cfg.method = 'realignfiducial';
  else
    cfg.method = 'interactive';
  end
end

% select the parameter that should be displayed
cfg.parameter = parameterselection(cfg.parameter, mri);
if iscell(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
end

switch cfg.method
  case 'realignfiducial'
    % do nothing
  case 'interactive'
    dat = getsubfield(mri, cfg.parameter);
    nas = [];
    lpa = [];
    rpa = [];
    x = 1:mri.dim(1);
    y = 1:mri.dim(2);
    z = 1:mri.dim(3);
    xc = round(mri.dim(1)/2);
    yc = round(mri.dim(2)/2);
    zc = round(mri.dim(3)/2);
    while(1) % break when 'q' is pressed
      fprintf('============================================================\n');
      fprintf('click with mouse button to reslice the display to a new position\n');
      fprintf('press n/l/r on keyboard to record the current position as fiducial location\n');
      fprintf('press q on keyboard to quit interactive mode\n');
      xc = round(xc);
      yc = round(yc);
      zc = round(zc);
      volplot(x, y, z, dat, [xc yc zc], cfg.clim);
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
      if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); else fprintf('nas = undefined\n'); end
      if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); else fprintf('lpa = undefined\n'); end
      if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); else fprintf('rpa = undefined\n'); end
    end

    cfg.fiducial.nas = nas;
    cfg.fiducial.lpa = lpa;
    cfg.fiducial.rpa = rpa;

  otherwise
    error('unsupported method');
end

% compute the homogenous transformation matrix describing the new coordinate system
vox2head  = headcoordinates(cfg.fiducial.nas, cfg.fiducial.lpa, cfg.fiducial.rpa);

if ~isfield(mri, 'transform')
  mri.transform = vox2head;
elseif all(all(mri.transform==eye(4)))
  mri.transform = vox2head;
else
  warning('removing old transformation matrix');
  scale          = eye(4);
  %FIXME check whether the following is ever necessary
  %origvox2head   = mri.transform;
  %scale(1:3,1:3) = diag(sqrt(sum(origvox2head(1:3,1:3).^2,2)));
  mri.transform  = scale*vox2head;
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: volumerealign.m,v 1.12 2009/07/31 13:43:36 jansch Exp $';

% remember the configuration
mri.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to show three orthogonal slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function volplot(x, y, z, dat, c, cscale);
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
crosshair([x(xi) z(zi)], 'color', 'yellow');

subplot(h2);
imagesc(y, z, squeeze(dat(xi,:,:))'); set(gca, 'ydir', 'normal')
axis equal; axis tight;
xlabel('j'); ylabel('k');
caxis([cmin cmax]);
crosshair([y(yi) z(zi)], 'color', 'yellow');

subplot(h3);
imagesc(x, y, squeeze(dat(:,:,zi))'); set(gca, 'ydir', 'normal')
axis equal; axis tight;
xlabel('i'); ylabel('j');
caxis([cmin cmax]);
crosshair([x(xi) y(yi)], 'color', 'yellow');

colormap gray
