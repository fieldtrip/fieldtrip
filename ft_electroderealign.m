function [norm] = ft_electroderealign(cfg)

% FT_ELECTRODEREALIGN rotates and translates electrode positions to
% template electrode positions or towards the head surface. It can
% either perform a rigid body transformation, in which only the
% coordinate system is changed, or it can apply additional deformations
% to the input electrodes.
%
% Use as
%   [elec] = ft_electroderealign(cfg)
% 
% This function requires the same electrode and channel positions.
% 
% Different methods for aligning the input electrodes to the subjects head
% are implemented, which are described in detail below:
%
% TEMPLATE - You can apply a spatial transformation/deformation that
% automatically minimizes the distance between the electrodes and the
% template or standard electrode set. The warping methods use a non-linear
% search to minimize the error between the input electrodes and
% corresponding template electrodes or between the input electrodes and a
% head surface.
%
% FIDUCIAL - You can apply a rigid body realignment based on three fiducial
% locations. Realigning using the fiducials only ensures that the fiducials
% (typically nose, left and right ear) are along the same axes in the input
% eectrode set as in the template electrode set.
%
% INTERACTIVE - You can display the skin surface together with the
% electrode position, and manually (using the graphical user interface)
% adjust the rotation, translation and scaling parameters, so that the
% electrodes correspond with the skin.
%
% MANUAL - You can display the skin surface and manually determine the
% electrode positions by clicking on the skin surface.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for aligning or placing the electrodes
%                        'template'        realign the electrodes to a template electrode set
%                        'fiducial'        realign using the NAS, LPA and RPA fiducials
%                        'interactive'     realign manually using a graphical user interface
%                        'manual'          manual positioning of the electrodes by clicking in a graphical user interface
%   cfg.warp          = string describing the spatial transformation for the template method
%                        'rigidbody'       apply a rigid-body warp (default)
%                        'globalrescale'   apply a rigid-body warp with global rescaling
%                        'traditional'     apply a rigid-body warp with individual axes rescaling
%                        'nonlin1'         apply a 1st order non-linear warp
%                        'nonlin2'         apply a 2nd order non-linear warp
%                        'nonlin3'         apply a 3rd order non-linear warp
%                        'nonlin4'         apply a 4th order non-linear warp
%                        'nonlin5'         apply a 5th order non-linear warp
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see  FT_CHANNELSELECTION for details
%   cfg.fiducial       = cell-array with the name of three fiducials used for
%                        realigning (default = {'nasion', 'lpa', 'rpa'})
%   cfg.casesensitive  = 'yes' or 'no', determines whether string comparisons
%                        between electrode labels are case sensitive (default = 'yes')
%   cfg.feedback       = 'yes' or 'no' (default = 'no')
%
% The EEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%
% If you want to align the electrodes to a single template electrode set
% or to multiple electrode sets (which will be averaged), you should
% specify the template electrode sets as
%   cfg.template       = single electrode set that serves as standard
% or
%   cfg.template{1..N} = list of electrode sets that are averaged into the standard
% The template electrode sets can be specified either as electrode
% structures (i.e. when they are already read in memory) or as electrode
% files.
%
% If you only want to realign using the fiducials, the template has to contain
% the three fiducials, e.g.
%   cfg.template.chanpos(1,:) = [110 0 0]  % location of the nose
%   cfg.template.chanpos(2,:) = [0  90 0]  % left ear
%   cfg.template.chanpos(3,:) = [0 -90 0]  % right ear
%   cfg.template.label    = {'nasion', 'lpa', 'rpa'}
%
% If you want to align existing electrodes to the head surface or position
% new electrodes on the head surface, you should specify the head surface as
%   cfg.headshape      = a filename containing headshape, a structure containing a
%                        single triangulated boundary, or a Nx3 matrix with surface
%                        points
%
% See also FT_READ_SENS, FT_VOLUMEREALIGN, FT_DATATYPE_SENS

% Copyright (C) 2005-2011, Robert Oostenveld
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
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

% text output
disp('Close the figure to output new sensor positions');

% set the defaults
if ~isfield(cfg, 'channel'),       cfg.channel = 'all';       end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'no';       end
if ~isfield(cfg, 'casesensitive'), cfg.casesensitive = 'yes'; end
if ~isfield(cfg, 'headshape'),     cfg.headshape = [];        end % for triangulated head surface, without labels
if ~isfield(cfg, 'template'),      cfg.template = [];         end % for electrodes or fiducials, always with labels
if ~isfield(cfg, 'warp'),          cfg.warp = 'rigidbody';    end
if ~isfield(cfg, 'label'),         cfg.label = 'off';        end % show labels

cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducials', 'fiducial'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'realignfiducial',  'fiducial'});
cfg = ft_checkconfig(cfg, 'forbidden', 'outline');
cfg = ft_checkconfig(cfg, 'renamedval',{'warp', 'homogenous', 'rigidbody'});

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% get the electrode definition that should be warped
if nargin==1
  try % try to get the description from the cfg
    elec = ft_fetch_sens(cfg);
  catch lasterr
    % start with an empty set of electrodes, this is useful for manual positioning
    elec = [];
    elec.pnt    = zeros(0,3);
    elec.label  = cell(0,1);
    elec.unit   = 'mm';
    warning(lasterr.message, lasterr.identifier);
  end
elseif nargin>1
  % the input electrodes were specified as second input argument
end
elec = ft_convert_units(elec); % ensure that the units are specified

% ensure up-to-date sensor descriptions 
elec = ft_datatype_sens(elec);

% ensure that channel and electrode positions are the same
assert(isequal(elec.elecpos,elec.chanpos),'This function requires same electrode and channel positions.'); 


usetemplate  = isfield(cfg, 'template')  && ~isempty(cfg.template);
useheadshape = isfield(cfg, 'headshape') && ~isempty(cfg.headshape);

if usetemplate
  % get the template electrode definitions
  if ~iscell(cfg.template)
    cfg.template = {cfg.template};
  end
  Ntemplate = length(cfg.template);
  for i=1:Ntemplate
    if isstruct(cfg.template{i})
      template(i) = cfg.template{i};
    else
      template(i) = ft_read_sens(cfg.template{i});
    end
  end
  
  clear tmp
  for i=1:Ntemplate
    tmp(i) = ft_convert_units(template(i), elec.unit); % ensure that the units are consistent with the electrodes
    tmp2(i) = ft_datatype_sens(tmp(i)); % ensure up-to-date sensor descriptions (Oct 2011)
  end
  template = tmp2;
end

if useheadshape
    % get the surface describing the head shape
    if isstruct(cfg.headshape) && isfield(cfg.headshape, 'hex')
        if isfield(cfg.headshape,'pos') && ~isfield(cfg.headshape,'pnt')
            cfg.headshape.pnt = cfg.headshape.pos;
        end
        headshape = mesh2edge(cfg.headshape);
    elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'tet')
        if isfield(cfg.headshape,'pos') && ~isfield(cfg.headshape,'pnt')
            cfg.headshape.pnt = cfg.headshape.pos;
        end
        headshape = mesh2edge(cfg.headshape);
    elseif isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
        % use the headshape surface specified in the configuration
        headshape = cfg.headshape;
    elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
        % use the headshape points specified in the configuration
        headshape.pnt = cfg.headshape;
    elseif ischar(cfg.headshape)
        % read the headshape from file
        headshape = ft_read_headshape(cfg.headshape);
    else
        error('cfg.headshape is not specified correctly')
    end
    if ~isfield(headshape, 'tri') && ~isfield(headshape,'poly')
        % generate a closed triangulation from the surface points
        headshape.pnt = unique(headshape.pnt, 'rows');
        headshape.tri = projecttri(headshape.pnt);
    end
    headshape = ft_convert_units(headshape, elec.unit); % ensure that the units are consistent with the electrodes
end

% remember the original electrode locations and labels
orig = elec;

% convert all labels to lower case for string comparisons
% this has to be done AFTER keeping the original labels and positions
if strcmp(cfg.casesensitive, 'no')
  for i=1:length(elec.label)
    elec.label{i} = lower(elec.label{i});
  end
  for j=1:length(template)
    for i=1:length(template(j).label)
      template(j).label{i} = lower(template(j).label{i});
    end
  end
end

if strcmp(cfg.feedback, 'yes')
  % create an empty figure, continued below...
  figure
  axis equal
  axis vis3d
  hold on
  xlabel('x')
  ylabel('y')
  zlabel('z')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.method, 'template') && usetemplate
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % determine electrode selection and overlapping subset for warping
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  for i=1:Ntemplate
    cfg.channel = ft_channelselection(cfg.channel, template(i).label);
  end
  
  % make subselection of electrodes
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label = elec.label(datsel);
  elec.chanpos   = elec.chanpos(datsel,:);
  for i=1:Ntemplate
    [cfgsel, datsel] = match_str(cfg.channel, template(i).label);
    template(i).label = template(i).label(datsel);
    template(i).chanpos   = template(i).chanpos(datsel,:);
  end
  
  % compute the average of the template electrode positions
  all = [];
  for i=1:Ntemplate
    all = cat(3, all, template(i).chanpos);
  end
  avg    = mean(all,3);
  stderr = std(all, [], 3);
  
  fprintf('warping electrodes to template... '); % the newline comes later
  [norm.chanpos, norm.m] = warp_optim(elec.chanpos, avg, cfg.warp);
  norm.label = elec.label;
  
  dpre  = mean(sqrt(sum((avg - elec.chanpos).^2, 2)));
  dpost = mean(sqrt(sum((avg - norm.chanpos).^2, 2)));
  fprintf('mean distance prior to warping %f, after warping %f\n', dpre, dpost);
  
  if strcmp(cfg.feedback, 'yes')
    % plot all electrodes before warping
    my_plot3(elec.chanpos, 'r.');
    my_plot3(elec.chanpos(1,:), 'r*');
    my_plot3(elec.chanpos(2,:), 'r*');
    my_plot3(elec.chanpos(3,:), 'r*');
    my_text3(elec.chanpos(1,:), elec.label{1}, 'color', 'r');
    my_text3(elec.chanpos(2,:), elec.label{2}, 'color', 'r');
    my_text3(elec.chanpos(3,:), elec.label{3}, 'color', 'r');
    
    % plot all electrodes after warping
    my_plot3(norm.chanpos, 'm.');
    my_plot3(norm.chanpos(1,:), 'm*');
    my_plot3(norm.chanpos(2,:), 'm*');
    my_plot3(norm.chanpos(3,:), 'm*');
    my_text3(norm.chanpos(1,:), norm.label{1}, 'color', 'm');
    my_text3(norm.chanpos(2,:), norm.label{2}, 'color', 'm');
    my_text3(norm.chanpos(3,:), norm.label{3}, 'color', 'm');
    
    % plot the template electrode locations
    my_plot3(avg,      'b.');
    my_plot3(avg(1,:), 'b*');
    my_plot3(avg(2,:), 'b*');
    my_plot3(avg(3,:), 'b*');
    my_text3(avg(1,:), norm.label{1}, 'color', 'b');
    my_text3(avg(2,:), norm.label{2}, 'color', 'b');
    my_text3(avg(3,:), norm.label{3}, 'color', 'b');
    
    % plot lines connecting the input/warped electrode locations with the template locations
    my_line3(elec.chanpos, avg, 'color', 'r');
    my_line3(norm.chanpos, avg, 'color', 'm');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'template') && useheadshape
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % determine electrode selection and overlapping subset for warping
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  
  % make subselection of electrodes
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label = elec.label(datsel);
  elec.chanpos   = elec.chanpos(datsel,:);
  
  fprintf('warping electrodes to head shape... '); % the newline comes later
  [norm.chanpos, norm.m] = warp_optim(elec.chanpos, headshape, cfg.warp);
  norm.label = elec.label;
  
  dpre  = warp_error([],     elec.chanpos, headshape, cfg.warp);
  dpost = warp_error(norm.m, elec.chanpos, headshape, cfg.warp);
  fprintf('mean distance prior to warping %f, after warping %f\n', dpre, dpost);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'fiducial')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % try to determine the fiducials automatically if not specified
  if ~isfield(cfg, 'fiducial')
    option1 = {'nasion' 'left' 'right'};
    option2 = {'nasion' 'lpa' 'rpa'};
    option3 = {'nz' 'lpa' 'rpa'};
    if length(match_str(elec.label, option1))==3
      cfg.fiducial = option1;
    elseif length(match_str(elec.label, option2))==3
      cfg.fiducial = option2;
    elseif length(match_str(elec.label, option3))==3
      cfg.fiducial = option3;
    else
      error('could not determine three fiducials, please specify cfg.fiducial')
    end
  end
  fprintf('using fiducials {''%s'', ''%s'', ''%s''}\n', cfg.fiducial{1}, cfg.fiducial{2}, cfg.fiducial{3});
  
  % determine electrode selection
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [cfgsel, datsel] = match_str(cfg.channel, elec.label);
  elec.label = elec.label(datsel);
  elec.chanpos   = elec.chanpos(datsel,:);
  
  if length(cfg.fiducial)~=3
    error('you must specify three fiducials');
  end
  
  % do case-insensitive search for fiducial locations
  nas_indx = match_str(lower(elec.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{3}));
  if length(nas_indx)~=1 || length(lpa_indx)~=1 || length(rpa_indx)~=1
    error('not all fiducials were found in the electrode set');
  end
  elec_nas = elec.chanpos(nas_indx,:);
  elec_lpa = elec.chanpos(lpa_indx,:);
  elec_rpa = elec.chanpos(rpa_indx,:);
  
  % FIXME change the flow in the remainder
  % if one or more template electrode sets are specified, then align to the average of those
  % if no template is specified, then align so that the fiducials are along the axis
  
  % find the matching fiducials in the template and average them
  templ_nas = [];
  templ_lpa = [];
  templ_rpa = [];
  for i=1:Ntemplate
    nas_indx = match_str(lower(template(i).label), lower(cfg.fiducial{1}));
    lpa_indx = match_str(lower(template(i).label), lower(cfg.fiducial{2}));
    rpa_indx = match_str(lower(template(i).label), lower(cfg.fiducial{3}));
    if length(nas_indx)~=1 || length(lpa_indx)~=1 || length(rpa_indx)~=1
      error(sprintf('not all fiducials were found in template %d', i));
    end
    templ_nas(end+1,:) = template(i).chanpos(nas_indx,:);
    templ_lpa(end+1,:) = template(i).chanpos(lpa_indx,:);
    templ_rpa(end+1,:) = template(i).chanpos(rpa_indx,:);
  end
  templ_nas = mean(templ_nas,1);
  templ_lpa = mean(templ_lpa,1);
  templ_rpa = mean(templ_rpa,1);
  
  % realign both to a common coordinate system
  elec2common  = headcoordinates(elec_nas, elec_lpa, elec_rpa);
  templ2common = headcoordinates(templ_nas, templ_lpa, templ_rpa);
  
  % compute the combined transform and realign the electrodes to the template
  norm       = [];
  norm.m     = elec2common * inv(templ2common);
  norm.chanpos   = warp_apply(norm.m, elec.chanpos, 'homogeneous');
  norm.label = elec.label;
  
  nas_indx = match_str(lower(elec.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(elec.label), lower(cfg.fiducial{3}));
  dpre  = mean(sqrt(sum((elec.chanpos([nas_indx lpa_indx rpa_indx],:) - [templ_nas; templ_lpa; templ_rpa]).^2, 2)));
  nas_indx = match_str(lower(norm.label), lower(cfg.fiducial{1}));
  lpa_indx = match_str(lower(norm.label), lower(cfg.fiducial{2}));
  rpa_indx = match_str(lower(norm.label), lower(cfg.fiducial{3}));
  dpost = mean(sqrt(sum((norm.chanpos([nas_indx lpa_indx rpa_indx],:) - [templ_nas; templ_lpa; templ_rpa]).^2, 2)));
  fprintf('mean distance between fiducials prior to realignment %f, after realignment %f\n', dpre, dpost);
  
  if strcmp(cfg.feedback, 'yes')
    % plot the first three electrodes before transformation
    my_plot3(elec.chanpos(1,:), 'r*');
    my_plot3(elec.chanpos(2,:), 'r*');
    my_plot3(elec.chanpos(3,:), 'r*');
    my_text3(elec.chanpos(1,:), elec.label{1}, 'color', 'r');
    my_text3(elec.chanpos(2,:), elec.label{2}, 'color', 'r');
    my_text3(elec.chanpos(3,:), elec.label{3}, 'color', 'r');
    
    % plot the template fiducials
    my_plot3(templ_nas, 'b*');
    my_plot3(templ_lpa, 'b*');
    my_plot3(templ_rpa, 'b*');
    my_text3(templ_nas, ' nas', 'color', 'b');
    my_text3(templ_lpa, ' lpa', 'color', 'b');
    my_text3(templ_rpa, ' rpa', 'color', 'b');
    
    % plot all electrodes after transformation
    my_plot3(norm.chanpos, 'm.');
    my_plot3(norm.chanpos(1,:), 'm*');
    my_plot3(norm.chanpos(2,:), 'm*');
    my_plot3(norm.chanpos(3,:), 'm*');
    my_text3(norm.chanpos(1,:), norm.label{1}, 'color', 'm');
    my_text3(norm.chanpos(2,:), norm.label{2}, 'color', 'm');
    my_text3(norm.chanpos(3,:), norm.label{3}, 'color', 'm');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'interactive')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % open a figure
  fig = figure;
  % add the data to the figure
  set(fig, 'CloseRequestFcn', @cb_close);
  setappdata(fig, 'elec', elec);
  setappdata(fig, 'transform', eye(4));
  if useheadshape
    setappdata(fig, 'headshape', headshape);
  end
  if usetemplate
    % FIXME interactive realigning to template electrodes is not yet supported
    % this requires a consistent handling of channel selection etc.
    setappdata(fig, 'template', template);
  end
  % add the GUI elements
  cb_creategui(gca);
  cb_redraw(gca);
  rotate3d on
  waitfor(fig);
  % get the data from the figure that was left behind as global variable
  global norm
  tmp = norm;
  clear global norm
  norm = tmp;
  clear tmp
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(cfg.method, 'manual')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % open a figure
  fig = figure;
  rotate3d on
  ft_plot_mesh(headshape, 'edgecolor', 'k')
  xyz = ft_select_point3d(headshape, 'multiple', true);
  orig.chanpos = xyz;
  for i=1:size(orig.chanpos,1)
    orig.label{i,1} = 'unknown';
  end
  
else
  error('unknown method');
end

% apply the spatial transformation to all electrodes, replace the
% electrode labels by their case-sensitive original values, and ensure the same
% electrode and channel positions
switch cfg.method
  case 'template'
    norm.chanpos   = warp_apply(norm.m, orig.chanpos, cfg.warp);
    norm.elecpos   = norm.chanpos;
  case {'fiducial' 'interactive'}
    norm.chanpos   = warp_apply(norm.m, orig.chanpos);
    norm.elecpos   = norm.chanpos;
  case 'manual'
    % the positions are already assigned in correspondence with the mesh
    norm = orig;
    norm.elecpos   = norm.chanpos;
  otherwise
    error('unknown method');
end

if isfield(orig, 'label')
  norm.label = orig.label;
end

norm = ft_datatype_sens(norm);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble history norm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some simple SUBFUNCTIONs that facilitate 3D plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = my_plot3(xyz, varargin)
h = plot3(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function h = my_text3(xyz, varargin)
h = text(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function my_line3(xyzB, xyzE, varargin)
for i=1:size(xyzB,1)
  line([xyzB(i,1) xyzE(i,1)], [xyzB(i,2) xyzE(i,2)], [xyzB(i,3) xyzE(i,3)], varargin{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to layout a moderately complex graphical user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = layoutgui(fig, geometry, position, style, string, value, tag, callback);
horipos  = geometry(1); % lower left corner of the GUI part in the figure
vertpos  = geometry(2); % lower left corner of the GUI part in the figure
width    = geometry(3); % width  of the GUI part in the figure
height   = geometry(4); % height of the GUI part in the figure
horidist = 0.05;
vertdist = 0.05;
options  = {'units', 'normalized', 'HorizontalAlignment', 'center'}; %  'VerticalAlignment', 'middle'
Nrow     = size(position,1);
h        = cell(Nrow,1);
for i=1:Nrow
  if isempty(position{i})
    continue;
  end
  position{i} = position{i} ./ sum(position{i});
  Ncol = size(position{i},2);
  ybeg = (Nrow-i  )/Nrow + vertdist/2;
  yend = (Nrow-i+1)/Nrow - vertdist/2;
  for j=1:Ncol
    xbeg    = sum(position{i}(1:(j-1))) + horidist/2;
    xend    = sum(position{i}(1:(j  ))) - horidist/2;
    pos(1) = xbeg*width  + horipos;
    pos(2) = ybeg*height + vertpos;
    pos(3) = (xend-xbeg)*width;
    pos(4) = (yend-ybeg)*height;
    h{i}{j} = uicontrol(fig, ...
      options{:}, ...
      'position', pos, ...
      'style',    style{i}{j}, ...
      'string',   string{i}{j}, ...
      'tag',      tag{i}{j}, ...
      'value',    value{i}{j}, ...
      'callback', callback{i}{j} ...
      );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(hObject, eventdata, handles);
% define the position of each GUI element
% position = {
%   [2 1 1 1]
%   [2 1 1 1]
%   [2 1 1 1]
%   [1]
%   [1]
%   [1]
%   [1]
%   [1 1]
%   };
% 
% % define the style of each GUI element
% style = {
%   {'text' 'edit' 'edit' 'edit'}
%   {'text' 'edit' 'edit' 'edit'}
%   {'text' 'edit' 'edit' 'edit'}
%   {'pushbutton'}
%   {'pushbutton'}
%   {'toggle'}
%   {'toggle'}
%   {'text' 'edit'}
%   };
% 
% % define the descriptive string of each GUI element
% string = {
%   {'rotate'    0 0 0}
%   {'translate' 0 0 0}
%   {'scale'     1 1 1}
%   {'redisplay'}
%   {'apply'}
%   {'toggle labels'}
%   {'toggle axes'}
%   {'alpha' 0.7}
%   };
% 
% % define the value of each GUI element
% value = {
%   {[] [] [] []}
%   {[] [] [] []}
%   {[] [] [] []}
%   {[]}
%   {[]}
%   {0}
%   {0}
%   {[] []}
%   };
% 
% % define a tag for each GUI element
% tag = {
%   {'' 'rx' 'ry' 'rz'}
%   {'' 'tx' 'ty' 'tz'}
%   {'' 'sx' 'sy' 'sz'}
%   {''}
%   {''}
%   {'toggle labels'}
%   {'toggle axes'}
%   {'' 'alpha'}
%   };
% 
% % define the callback function of each GUI element
% callback = {
%   {[] @cb_redraw @cb_redraw @cb_redraw}
%   {[] @cb_redraw @cb_redraw @cb_redraw}
%   {[] @cb_redraw @cb_redraw @cb_redraw}
%   {@cb_redraw}
%   {@cb_apply}
%   {@cb_redraw}
%   {@cb_redraw}
%   {[] @cb_redraw}
%   };
% 
% fig = get(hObject, 'parent');
% layoutgui(fig, [0.7 0.05 0.25 0.50], position, style, string, value, tag, callback);

fig = get(hObject, 'parent');
% constants
CONTROL_WIDTH  = 0.05;
CONTROL_HEIGHT = 0.06; 
CONTROL_HOFFSET = 0.7;
CONTROL_VOFFSET = 0.5;
% rotateui
uicontrol('tag', 'rotateui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'rotate', 'callback', [])
uicontrol('tag', 'rx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ry', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'rz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'rotateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'rx',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'ry',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'rz',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);

% scaleui
uicontrol('tag', 'scaleui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'scale', 'callback', [])
uicontrol('tag', 'sx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sy', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'scaleui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sx',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sy',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sz',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ty', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'tz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tx',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ty',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tz',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);

% control buttons
uicontrol('tag', 'redisplaybtn', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'redisplay', 'value', [], 'callback', @cb_redraw);
uicontrol('tag', 'applybtn', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'apply',     'value', [], 'callback', @cb_apply);
uicontrol('tag', 'toggle labels', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle label', 'value', 0, 'callback', @cb_redraw);
uicontrol('tag', 'toggle axes', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle axes', 'value', 0, 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'redisplaybtn',  'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-3*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'applybtn',      'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle labels', 'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle axes',   'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);

% alpha ui (somehow not implemented, facealpha is fixed at 0.7
uicontrol('tag', 'alphaui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'alpha', 'value', [], 'callback', []);
uicontrol('tag', 'alpha',   'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0.7',   'value', [], 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'alphaui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 3*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'alpha',   'BackgroundColor', [0.8 0.8 0.8], 'width', 3*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(hObject, eventdata, handles);
fig = get(hObject, 'parent');
headshape = getappdata(fig, 'headshape');
bnd.pnt = headshape.pnt; %ft_plot_mesh wants headshape in bnd fields
if isfield(headshape,'tri')
bnd.tri = headshape.tri;
elseif isfield(headshape, 'poly')
bnd.poly = headshape.poly;
end
elec = getappdata(fig, 'elec');
template = getappdata(fig, 'template');
% get the transformation details
rx = str2num(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2num(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2num(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2num(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2num(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2num(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2num(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2num(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2num(get(findobj(fig, 'tag', 'sz'), 'string'));
R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
elec = ft_transform_sens(H, elec);
axis vis3d; cla
xlabel('x')
ylabel('y')
zlabel('z')

if ~isempty(template)
  if size(template.chanpos, 2)==2
    hs = plot(template.chanpos(:,1), template.chanpos(:,2), 'b.', 'MarkerSize', 20);
  else
    hs = plot3(template.chanpos(:,1), template.chanpos(:,2), template.chanpos(:,3), 'b.', 'MarkerSize', 20);
  end
end

if ~isempty(headshape)
  % plot the faces of the 2D or 3D triangulation
  skin   = [255 213 119]/255;
  brain  = [202 100 100]/255;
  cortex = [255 213 119]/255;
  ft_plot_mesh(bnd,'facecolor', skin,'EdgeColor','none','facealpha',0.7);
  lighting gouraud
  material shiny
  camlight
end


if isfield(elec, 'fid') && ~isempty(elec.fid.chanpos)
  ft_plot_sens(elec.fid,'style', 'r*');
end
if get(findobj(fig, 'tag', 'toggle axes'), 'value')
  axis on
  grid on
else
  axis off
  grid on
end

if get(findobj(fig, 'tag', 'toggle labels'), 'value')
  cfg.label = 'on';
else
  cfg.label = 'off';
end
hold on
ft_plot_sens(elec,'label',cfg.label);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_apply(hObject, eventdata, handles);
fig = get(hObject, 'parent');
elec      = getappdata(fig, 'elec');
transform = getappdata(fig, 'transform');
% get the transformation details
rx = str2num(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2num(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2num(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2num(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2num(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2num(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2num(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2num(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2num(get(findobj(fig, 'tag', 'sz'), 'string'));
R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
elec = ft_transform_headshape(H, elec);
transform = H * transform;
set(findobj(fig, 'tag', 'rx'), 'string', 0);
set(findobj(fig, 'tag', 'ry'), 'string', 0);
set(findobj(fig, 'tag', 'rz'), 'string', 0);
set(findobj(fig, 'tag', 'tx'), 'string', 0);
set(findobj(fig, 'tag', 'ty'), 'string', 0);
set(findobj(fig, 'tag', 'tz'), 'string', 0);
set(findobj(fig, 'tag', 'sx'), 'string', 1);
set(findobj(fig, 'tag', 'sy'), 'string', 1);
set(findobj(fig, 'tag', 'sz'), 'string', 1);
setappdata(fig, 'elec', elec);
setappdata(fig, 'transform', transform);
cb_redraw(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_close(hObject, eventdata, handles);
% make the current transformation permanent and subsequently allow deleting the figure
cb_apply(gca);
% get the updated electrode from the figure
fig    = hObject;
% hmmm, this is ugly
global norm
norm   = getappdata(fig, 'elec');
norm.m = getappdata(fig, 'transform');
set(fig, 'CloseRequestFcn', @delete);
delete(fig);

