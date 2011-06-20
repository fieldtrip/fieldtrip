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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_READ_MRI, FT_ELECTRODEREALIGN

% Undocumented option
%   cfg.coordsys, works for interactive and fiducial

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

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial', 'fiducial'});
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'required', 'method');

% set the defaults
if ~isfield(cfg, 'fiducial'),  cfg.fiducial = [];         end
if ~isfield(cfg, 'landmark'),  cfg.landmark = [];         end
if ~isfield(cfg, 'parameter'), cfg.parameter = 'anatomy'; end
if ~isfield(cfg, 'clim'),      cfg.clim      = [];        end
if ~isfield(cfg, 'inputfile'), cfg.inputfile = [];        end
if ~isfield(cfg, 'outputfile'),cfg.outputfile = [];       end
if ~isfield(cfg, 'coordsys') && (strcmp(cfg.method, 'interactive') || strcmp(cfg.method, 'fiducial'))
  cfg.coordsys = 'ctf';
end

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
elseif strcmp(cfg.method, 'fiducial')
  basedonfid = 1;
  basedonmrk = 0;
elseif strcmp(cfg.method, 'landmark')
  basedonfid = 0;
  basedonmrk = 1;
end

% select the parameter that should be displayed
cfg.parameter = parameterselection(cfg.parameter, mri);
if iscell(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
end

h1 = subplot('position',[0.02 0.55 0.44 0.44]);%subplot(2,2,1);
h2 = subplot('position',[0.52 0.55 0.44 0.44]);%subplot(2,2,2);
h3 = subplot('position',[0.02 0.05 0.44 0.44]);%subplot(2,2,3);
handles = {h1 h2 h3};
switch cfg.method
  case 'fiducial'
    % do nothing
    
  case 'landmark'
    % do nothing
    
  case 'interactive'
    showcrosshair = true;
    showmarkers   = true;
    dat = mri.(cfg.parameter);
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
    
    updatepanel = [1 2 3];
    pnt         = zeros(0,3);
    markerpos   = zeros(0,3);
    markerlabel = {};
    markercolor = {};
    while true % break when 'q' is pressed
      fprintf('click with mouse button to reslice the display to a new position\n');
      fprintf('press n/l/r on keyboard to record the current position as fiducial location\n');
      fprintf('press a/p/z on keyboard to record the current position as anatomical landmark\n');
      fprintf('press the arrow keys on the keyboard to increment or decrement the slice number by one\n');
      fprintf('press c or C on the keyboard to show or hide the crosshair\n');
      fprintf('press q on keyboard to quit interactive mode\n');
      
      xc = round(xc); xc = max(1,xc); xc = min(mri.dim(1),xc);
      yc = round(yc); yc = max(1,yc); yc = min(mri.dim(2),yc);
      zc = round(zc); zc = max(1,zc); zc = min(mri.dim(3),zc);
      markers = {markerpos markerlabel markercolor};
      [h1, h2, h3] = volplot(x, y, z, dat, [xc yc zc], cfg.clim, showcrosshair, updatepanel, handles, showmarkers, markers);
      drawnow;
      try, [d1, d2, key] = ginput(1); catch, key='q'; end
      switch key
      case 113 % 'q'
        break;
      case 108 % 'l'
        lpa = [xc yc zc];
      case 114 % 'r'
        rpa = [xc yc zc];
      case 110 % 'n'
        nas = [xc yc zc];
      case 97  % 'a'
        antcomm = [xc yc zc];
      case 112 % 'p'
        pstcomm = [xc yc zc];
      case 122 % 'z'
        xzpoint = [xc yc zc];
      case 99  % 'c'
        showcrosshair = true;
      case 67  % 'C'
        showcrosshair = false;
      case 1 % left mouse click
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
          otherwise
            continue;
        end
        switch l2,
          case 'i'
            xc = d2;
          case 'j'
            yc = d2;
          case 'k'
            zc = d2;
          otherwise
            continue;
        end
        if l1=='i' && l2=='j'
          updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          updatepanel = [3 1 2];
        end
      case 3 % right mouse click  
        % add point to a list
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
        pnt = [pnt; xc yc zc];
        if l1=='i' && l2=='j'
          updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          updatepanel = [3 1 2];
        end
      case 2 % middle mouse click
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        
        % remove the previous point
        if size(pnt,1)>0
          pnt(end,:) = [];
        end
        
        if l1=='i' && l2=='j'
          updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          updatepanel = [3 1 2];
        end
      case 28 % arrow left
        % update the coordinates
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        if l1=='i' && l2=='j'
          xc = xc-1; updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          xc = xc-1; updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          yc = yc-1; updatepanel = [3 1 2];
        end
      case 30 % arrow up
        % update the coordinates
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        if l1=='i' && l2=='j'
          yc = yc+1; updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          zc = zc+1; updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          zc = zc+1; updatepanel = [3 1 2];
        end
      case 29 % arrow right
        % update the coordinates
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        if l1=='i' && l2=='j'
          xc = xc+1; updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          xc = xc+1; updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          yc = yc+1; updatepanel = [3 1 2];
        end
      case 31 % arrow down
        % update the coordinates
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        if l1=='i' && l2=='j'
          yc = yc-1; updatepanel = [1 2 3];
        elseif l1=='i' && l2=='k'
          zc = zc-1; updatepanel = [2 3 1];
        elseif l1=='j' && l2=='k'
          zc = zc-1; updatepanel = [3 1 2];
        end
      otherwise
        % do nothing
      end
      fprintf('============================================================\n');
      
      if all(round([xc yc zc])<=mri.dim)
        str = sprintf('voxel %d, indices [%d %d %d]', sub2ind(mri.dim(1:3), round(xc), round(yc), round(zc)), round([xc yc zc]));
            
        if isfield(mri, 'coordsys') && isfield(mri, 'unit')
          str = sprintf('%s, %s coordinates [%.1f %.1f %.1f] %s', str, mri.coordsys, warp_apply(mri.transform, [xc yc zc]), mri.unit);
        elseif ~isfield(mri, 'coordsys') && isfield(mri, 'unit')
          str = sprintf('%s, location [%.1f %.1f %.1f] %s', str, warp_apply(mri.transform, [xc yc zc]), mri.unit);
        elseif isfield(mri, 'coordsys') && ~isfield(mri, 'unit')
          str = sprintf('%s, %s coordinates [%.1f %.1f %.1f]', str, mri.coordsys, warp_apply(mri.transform, [xc yc zc]));
        elseif ~isfield(mri, 'coordsys') && ~isfield(mri, 'unis')
          str = sprintf('%s, location [%.1f %.1f %.1f]', str, warp_apply(mri.transform, [xc yc zc]));
        end
        fprintf('%s\n', str);
        % fprintf('cur_voxel = [%f %f %f], cur_head = [%f %f %f]\n', [xc yc zc], warp_apply(mri.transform, [xc yc zc]));
      end
      
      markerpos   = zeros(0,3);
      markerlabel = {};
      markercolor = {};
      if ~isempty(nas),
        fprintf('nas_voxel = [%f %f %f], nas_head = [%f %f %f]\n', nas, warp_apply(mri.transform, nas));
        markerpos   = [markerpos; nas];
        markerlabel = [markerlabel; {'nas'}];
        markercolor = [markercolor; {'b'}];
      end
      if ~isempty(lpa),
        fprintf('lpa_voxel = [%f %f %f], lpa_head = [%f %f %f]\n', lpa, warp_apply(mri.transform, lpa));
        markerpos   = [markerpos; lpa];
        markerlabel = [markerlabel; {'lpa'}];
        markercolor = [markercolor; {'g'}];
      end
      if ~isempty(rpa),
        fprintf('rpa_voxel = [%f %f %f], rpa_head = [%f %f %f]\n', rpa, warp_apply(mri.transform, rpa));
        markerpos   = [markerpos; rpa];
        markerlabel = [markerlabel; {'rpa'}];
        markercolor = [markercolor; {'r'}];
      end
      if ~isempty(antcomm),
        fprintf('antcomm_voxel = [%f %f %f], antcomm_head = [%f %f %f]\n', antcomm, warp_apply(mri.transform, antcomm));
        markerpos   = [markerpos; antcomm];
        markerlabel = [markerlabel; {'antcomm'}];
        markercolor = [markercolor; {'b'}];
      end
      if ~isempty(pstcomm),
        fprintf('pstcomm_voxel = [%f %f %f], pstcomm_head = [%f %f %f]\n', pstcomm, warp_apply(mri.transform, pstcomm));
        markerpos   = [markerpos; pstcomm];
        markerlabel = [markerlabel; {'pstcomm'}];
        markercolor = [markercolor; {'g'}];
      end
      if ~isempty(xzpoint),
        fprintf('xzpoint_voxel = [%f %f %f], xzpoint_head = [%f %f %f]\n', xzpoint, warp_apply(mri.transform, xzpoint));
        markerpos   = [markerpos; xzpoint];
        markerlabel = [markerlabel; {'xzpoint'}];
        markercolor = [markercolor; {'r'}];
      end
      if ~isempty(pnt)
        fprintf('%f extra points selected\n', size(pnt,1));
        markerpos   = [markerpos; pnt];
        markerlabel = [markerlabel; repmat({''}, size(pnt,1), 1)];
        markercolor = [markercolor; repmat({'m'}, size(pnt,1), 1)];
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
  [realign, coordsys] = headcoordinates(nas_head, lpa_head, rpa_head, cfg.coordsys);
  
elseif basedonmrk
  % the fiducial locations are now specified in voxels, convert them to head
  % coordinates according to the existing transform matrix
  ac     = warp_apply(mri.transform, cfg.landmark.ac);
  pc     = warp_apply(mri.transform, cfg.landmark.pc);
  xzpoint= warp_apply(mri.transform, cfg.landmark.xzpoint);
  
  % compute the homogenous transformation matrix describing the new coordinate system
  [realign, coordsys] = headcoordinates(ac, pc, xzpoint, 'spm');
  
else
  realign = [];
  
end

if ~isempty(realign)
  % combine the additional transformation with the original one
  mri.transformorig = mri.transform;
  mri.transform = realign * mri.transform;
  mri.coordsys  = coordsys;
else
  warning('no coordinate system reallignment has been done');
end

if exist('pnt', 'var')
  mri.pnt = pnt;
end

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
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration
mri.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'mri', mri); % use the variable name "data" in the output file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to show three orthogonal slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h1, h2, h3] = volplot(x, y, z, dat, c, cscale, showcrosshair, updatepanel, handles, showmarkers, markers)

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

if nargin<8
  updatepanel = [1 2 3];
end

if nargin<9
  h1 = [];
  h2 = [];
  h3 = [];
else
  h1 = handles{1};
  h2 = handles{2};
  h3 = handles{3};
end

if showmarkers
  markerpos = round(markers{1});
  markercolor = markers{3};
  sel1 = find(markerpos(:,2)==repmat(c(2),size(markerpos,1),1));
  sel2 = find(markerpos(:,1)==repmat(c(1),size(markerpos,1),1));
  sel3 = find(markerpos(:,3)==repmat(c(3),size(markerpos,1),1));
end

for k = 1:numel(updatepanel)

  update = updatepanel(k);
  if update==1
    subplot(h1);
    imagesc(x, z, squeeze(dat(:,yi,:))'); set(gca, 'ydir', 'normal')
    axis equal; axis tight;
    xlabel('i'); ylabel('k');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([x(xi) z(zi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel1)>0
      hold on;
      for kk = 1:numel(sel1)
        plot(markerpos(sel1(kk),1), markerpos(sel1(kk),3), 'marker', '.', 'color', markercolor{sel1(kk)});
      end
      hold off;
    end
  end
  
  if update==2
    subplot(h2);
    imagesc(y, z, squeeze(dat(xi,:,:))'); set(gca, 'ydir', 'normal')
    axis equal; axis tight;
    xlabel('j'); ylabel('k');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([y(yi) z(zi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel2)>0
      hold on;
      for kk = 1:numel(sel2)
        plot(markerpos(sel2(kk),2), markerpos(sel2(kk),3), 'marker', '.', 'color', markercolor{sel2(kk)});
      end
      hold off;
    end  
  end
  
  if update==3
    subplot(h3); 
    imagesc(x, y, squeeze(dat(:,:,zi))'); set(gca, 'ydir', 'normal')
    axis equal; axis tight;
    xlabel('i'); ylabel('j');
    caxis([cmin cmax]);
    if showcrosshair
      crosshair([x(xi) y(yi)], 'color', 'yellow');
    end
    if showmarkers && numel(sel3)>0
      hold on;
      for kk = 1:numel(sel3)
        plot(markerpos(sel3(kk),1), markerpos(sel3(kk),2), 'marker', '.', 'color', markercolor{sel3(kk)});
      end
      hold off;
    end
  end

end

colormap gray

h = gca;
