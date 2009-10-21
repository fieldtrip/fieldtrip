function [lay] = prepare_layout(cfg, data);

% PREPARE_LAYOUT creates a 2-D layout of the channel locations. This layout
% is required for plotting the topographical distribution of the potential
% or field distribution, or for plotting timecourses in a topographical
% arrangement.
%
% Use as
%   lay = prepare_layout(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.lay file, it can be created based on 3-D electrode or
% gradiometer positions in the configuration or in the data, or it can be
% created based on the specification of an electrode of gradiometer file.
%
% You can specify either one of the following configuration options
%   cfg.layout      filename containg the layout
%   cfg.rotate      number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  string, 2D projection method can be 'stereographic', 'orthographic', 'polar', 'gnomic' or 'inverse' (default = 'polar')
%   cfg.elec        structure with electrode positions, or
%   cfg.elecfile    filename containing electrode positions
%   cfg.grad        structure with gradiometer definition, or
%   cfg.gradfile    filename containing gradiometer definition
%   cfg.output      filename to which the layout will be written (default = [])
%   cfg.montage     'no' or a montage structure (default = 'no')
%   cfg.image       filename, use an image to construct a layout (e.g. usefull for ECoG grids)
%   cfg.bw          if an image is used and bw = 1 transforms the image in black and white (default = 0, do not transform)
%
% Alternatively the layout can be constructed from either
%   data.elec     structure with electrode positions
%   data.grad     structure with gradiometer definition
%
% Alternatively you can specify the following layouts which will be
% generated for all channels present in the data. Note that these layouts
% are suitable for multiplotting, but not for topoplotting.
%   cfg.layout = 'ordered'  will give you a NxN ordered layout
%   cfg.layout = 'vertical' will give you a Nx1 ordered layout
%   cfg.layout = 'butterfly'  will give you a layout with all channels on top of each other
%
% See also layoutplot, topoplotER, topoplotTFR, multiplotER, multiplotTFR

% TODO switch to using planarchannelset function

% undocumented and non-recommended option (for SPM only)
%   cfg.style       string, '2d' or '3d' (default = '2d')

% Copyright (C) 2007-2009, Robert Oostenveld
%
% $Log: prepare_layout.m,v $
% Revision 1.41  2009/10/14 15:26:36  roboos
% There was a typo in the cfg.projection default, which was ortographic
% (without "h").  The consequence was that the actual default being
% used is polar (because the underlying function would not recognise
% the typo and switch to its own polar default). I have changed the
% default to polar, consistent with its actual behaviour sofar. Polar
% is also the most robust for electrodes along the edge.
%
% Revision 1.40  2009/08/05 08:22:09  roboos
% better detection of empty/absent input data
%
% Revision 1.39  2009/08/05 06:32:41  roboos
% fixed layout generation for ordered and vertical when no data was given as input, labels fully depend on cfg.channel, not on data.label
%
% Revision 1.38  2009/08/04 13:57:00  roboos
% make tight vertical and orderer layout in case cfg.channel is specified
% allow skipping the COMNT and SCALE positions through the cfg
%
% Revision 1.37  2009/06/30 07:08:55  roboos
% removed debug keyboard statement
%
% Revision 1.36  2009/06/17 14:03:41  roboos
% consistent handling of ginput in case figure is closed
%
% Revision 1.35  2009/06/05 15:30:05  crimic
% updated help
%
% Revision 1.34  2009/06/05 15:28:03  crimic
% updated cfg
%
% Revision 1.33  2009/06/05 15:26:32  crimic
% minor change
%
% Revision 1.32  2009/05/18 15:59:44  roboos
% optinal plotting of RGB as greyscale image
%
% Revision 1.31  2009/03/30 17:55:17  vlalit
% Changed prepare_layout and headmodelplot to use channelposition. Changed the color
%  of sensor markers in headmodelplot to green for consistency with SPM plots.
%
% Revision 1.30  2009/02/26 10:37:10  roboos
% fixed bug, thanks to Yoni
%
% Revision 1.29  2009/02/25 17:44:27  roboos
% added vertical and butterfly layouts
%
% Revision 1.28  2009/02/04 16:58:31  roboos
% fixed backward compatibility bug in related to isfield, thanks to Irina
%
% Revision 1.27  2009/01/21 10:02:58  roboos
% explicit handling of config->struct in case configuration tracking is used
% changed default handling, don't use try-catch any more (thanks to Ivar)
%
% Revision 1.26  2009/01/20 14:27:24  jansch
% fixed bug in the case that the cfg.layout is a config-object leading to the
% generation of a 'default' CTF151 layout. this was problematic in the
% interactive mode of the plotting routines, toggling back and forth
%
% Revision 1.25  2008/12/24 13:34:19  roboos
% fixed bug in  hte closing of masking polygons
%
% Revision 1.24  2008/12/05 10:58:17  roboos
% fixed bug due to non-aligned vectors Rem and Lbl in subfunction readlay (thanks to pawel)
%
% Revision 1.23  2008/10/22 07:25:01  roboos
% undo the balancing of ctf meg channels prior to removing the unused coils
%
% Revision 1.22  2008/10/21 20:36:21  roboos
% (in grad2lay for ctf) prune the gradiometers to only contain meg
% channels and corresponding meg coils: this will not work for balanced
% gradiometers.
%
% Revision 1.21  2008/09/25 15:22:46  roboos
% fixed bug due to overlapping use of cfg.style, thanks to Tineke
%
% Revision 1.20  2008/09/25 12:55:52  roboos
% prevent x and y from being scaled in case of 3d layout, also don't add width, height, mask and outline for 3d
%
% Revision 1.19  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.18  2008/09/22 12:54:22  roboos
% added default handling for construction of outline and mask and for scaling the electrode positions to fith within unit circle
%
% Revision 1.17  2008/07/16 09:27:36  roboos
% added support for creating a layout based on a figure, including mask and outline for topoplot
%
% Revision 1.16  2008/06/25 06:36:06  roboos
% added option cfg.montage, computes average location for bipolar channels
%
% Revision 1.15  2008/05/14 19:17:24  roboos
% fixed bug for vladimirs neuromag306 example, removed spm specific documentation
%
% Revision 1.14  2008/05/14 13:53:30  roboos
% rotate only if needed
%
% Revision 1.13  2008/05/13 20:19:39  roboos
% changed senstype eeg into electrode
%
% Revision 1.12  2008/05/13 09:54:07  roboos
% added option cfg.style=2d|3d, used by SPM8
%
% Revision 1.11  2008/05/06 13:16:56  roboos
% added option for writing *.lay files
% merged bti and ctf code
%
% Revision 1.10  2008/04/25 12:29:52  roboos
% slight improvement for ordered layout
%
% Revision 1.9  2008/03/05 10:46:36  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.8  2007/12/12 09:59:10  roboos
% use cfg.feedback instead of global fb variable
%
% Revision 1.7  2007/11/05 09:43:31  roboos
% only whitespace
%
% Revision 1.6  2007/05/06 09:06:37  roboos
% implemented layout=ordered, for multiplotting only
%
% Revision 1.5  2007/03/21 15:52:30  roboos
% included the cfg.layout=lay case in the if-ladder
%
% Revision 1.4  2007/03/21 14:17:28  chrhes
% added a check to detect the (unlikely) case where cfg.layout already contains
% a valid layout (lay) structure, which is then returned as is; added a few
% comments to code; updated documentation.
%
% Revision 1.3  2007/03/20 10:41:30  roboos
% added options cfg.rotate and cfg.projection
% changed the default projection method from stereographic into orthographic
% by default rotate MEG electrode positions 90 degrees around the z-axis
% changed the default rotation for some MEG systems (now it is more explicit in
% the code)
%
% Revision 1.2  2007/03/20 09:23:23  chrhes
% small change to subfunction grad2lay which allows for MEG channel labels that
% do not have a space after the "MEG" token in the cases of neuromag122 and
% neuromag306 data
%
% Revision 1.1  2007/03/14 08:44:29  roboos
% new function that replaces private/createlayout, this new function can be used
% by end-users added support for mat files containing a lay variable, made some
% changes to the lay structure
%

% Undocumented option:
% cfg.layout can contain a lay structure which is simply returned as is

fieldtripdefs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1) || (nargin>2), error('incorrect number of input arguments'); end;
if (nargin<2), data = []; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default configuration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg, 'rotate'),     cfg.rotate = [];                end  % [] => rotation is determined based on the type of sensors
if ~isfield(cfg, 'style'),      cfg.style = '2d';               end
if ~isfield(cfg, 'projection'), cfg.projection = 'polar';       end
if ~isfield(cfg, 'layout'),     cfg.layout = [];                end
if ~isfield(cfg, 'grad'),       cfg.grad = [];                  end
if ~isfield(cfg, 'elec'),       cfg.elec = [];                  end
if ~isfield(cfg, 'gradfile'),   cfg.gradfile = [];              end
if ~isfield(cfg, 'elecfile'),   cfg.elecfile = [];              end
if ~isfield(cfg, 'output'),     cfg.output = [];                end
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'no';            end
if ~isfield(cfg, 'montage'),    cfg.montage = 'no';             end
if ~isfield(cfg, 'image'),      cfg.image = [];                 end
if ~isfield(cfg, 'bw'),         cfg.bw = 0;                     end
if ~isfield(cfg, 'channel'),    cfg.channel = 'all';            end 
if ~isfield(cfg, 'skipscale'),  cfg.skipscale = 'no';           end 
if ~isfield(cfg, 'skipcomnt'),  cfg.skipcomnt = 'no';           end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to generate the layout structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


skipscale = strcmp(cfg.skipscale, 'yes'); % in general a scale is desired
skipcomnt = strcmp(cfg.skipcomnt, 'yes'); % in general a comment desired

if isa(cfg.layout, 'config')
  % convert the nested config-object back into a normal structure
  cfg.layout = struct(cfg.layout);
end

% check whether cfg.layout already contains a valid layout structure (this can
% happen when higher level plotting functions are called with cfg.layout set to
% a lay structure)
if isstruct(cfg.layout) && isfield(cfg.layout, 'pos') && isfield(cfg.layout, 'label') && isfield(cfg.layout, 'width') && isfield(cfg.layout, 'height')
  lay = cfg.layout;

elseif isstruct(cfg.layout) && isfield(cfg.layout, 'pos') && isfield(cfg.layout, 'label') && (~isfield(cfg.layout, 'width') || ~isfield(cfg.layout, 'height'))
  lay = cfg.layout;
  % add width and height for multiplotting
  d = dist(lay.pos');
  nchans = length(lay.label);
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  lay.width  = ones(nchans,1) * mindist * 0.8;
  lay.height = ones(nchans,1) * mindist * 0.6;

elseif isequal(cfg.layout, 'butterfly')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    cfg.channel = channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    lay.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    lay.label = cfg.channel;
  end
  lay.pos     = zeros(nchan,2);  % centered at (0,0)
  lay.width   = ones(nchan,1) * 1.0;
  lay.height  = ones(nchan,1) * 1.0;
  lay.mask    = {};
  lay.outline = {};
  skipscale = true; % a scale is not desired
  skipcomnt = true; % a comment is initially not desired, or at least requires more thought

elseif isequal(cfg.layout, 'vertical')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    cfg.channel = channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    lay.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    lay.label = cfg.channel;
  end
  for i=1:(nchan+2)
    x = 0.5;
    y = 1-i/(nchan+1+2);
    lay.pos   (i,:) = [x y];
    lay.width (i,1) = 0.9;
    lay.height(i,1) = 0.9 * 1/(nchan+1+2);
    if i==(nchan+1)
      lay.label{i}   = 'SCALE';
    elseif i==(nchan+2)
      lay.label{i}   = 'COMNT';
    end
  end
  lay.mask    = {};
  lay.outline = {};

elseif isequal(cfg.layout, 'ordered')
  if nargin>1 && ~isempty(data)
    % look at the data to determine the overlapping channels
    cfg.channel = channelselection(cfg.channel, data.label);
    chanindx    = match_str(data.label, cfg.channel);
    nchan       = length(data.label(chanindx));
    lay.label   = data.label(chanindx);
  else
    nchan     = length(cfg.channel);
    lay.label = cfg.channel;
  end
  ncol = ceil(sqrt(nchan))+1;
  nrow = ceil(sqrt(nchan))+1;
  k = 0;
  for i=1:nrow
    for j=1:ncol
      k = k+1;
      if k<=nchan
        x = (j-1)/ncol;
        y = (nrow-i-1)/nrow;
        lay.pos(k,:) = [x y];
        lay.width(k,1)  = 0.8 * 1/ncol;
        lay.height(k,1) = 0.8 * 1/nrow;
      end
    end
  end

  lay.label{end+1}  = 'SCALE';
  lay.width(end+1)  = 0.8 * 1/ncol;
  lay.height(end+1) = 0.8 * 1/nrow;
  x = (ncol-2)/ncol;
  y = 0/nrow;
  lay.pos(end+1,:) = [x y];

  lay.label{end+1}  = 'COMNT';
  lay.width(end+1)  = 0.8 * 1/ncol;
  lay.height(end+1) = 0.8 * 1/nrow;
  x = (ncol-1)/ncol;
  y = 0/nrow;
  lay.pos(end+1,:) = [x y];

  % try to generate layout from other configuration options
elseif ischar(cfg.layout) && filetype(cfg.layout, 'matlab')
  fprintf('reading layout from file %s\n', cfg.layout);
  load(cfg.layout, 'lay');

elseif ischar(cfg.layout) && filetype(cfg.layout, 'layout')
  fprintf('reading layout from file %s\n', cfg.layout);
  lay = readlay(cfg.layout);

elseif ischar(cfg.layout) && ~filetype(cfg.layout, 'layout')
  % assume that cfg.layout is an electrode file
  fprintf('creating layout from electrode file %s\n', cfg.layout);
  lay = sens2lay(read_sens(cfg.layout), cfg.rotate, cfg.projection, cfg.style);

elseif ischar(cfg.elecfile)
  fprintf('creating layout from electrode file %s\n', cfg.elecfile);
  lay = sens2lay(read_sens(cfg.elecfile), cfg.rotate, cfg.projection, cfg.style);

elseif ~isempty(cfg.elec) && isstruct(cfg.elec)
  fprintf('creating layout from cfg.elec\n');
  lay = sens2lay(cfg.elec, cfg.rotate, cfg.projection, cfg.style);

elseif isfield(data, 'elec') && isstruct(data.elec)
  fprintf('creating layout from data.elec\n');
  lay = sens2lay(data.elec, cfg.rotate, cfg.projection, cfg.style);

elseif ischar(cfg.gradfile)
  fprintf('creating layout from gradiometer file %s\n', cfg.gradfile);
  lay = sens2lay(read_sens(cfg.gradfile), cfg.rotate, cfg.projection, cfg.style);

elseif ~isempty(cfg.grad) && isstruct(cfg.grad)
  fprintf('creating layout from cfg.grad\n');
  lay = sens2lay(cfg.grad, cfg.rotate, cfg.projection, cfg.style);

elseif isfield(data, 'grad') && isstruct(data.grad)
  fprintf('creating layout from data.grad\n');
  lay = sens2lay(data.grad, cfg.rotate, cfg.projection, cfg.style);

elseif ~isempty(cfg.image) && isempty(cfg.layout)
  fprintf('reading background image from %s\n', cfg.image);
  img = imread(cfg.image);
  img = flipdim(img, 1); % in combination with "axis xy"

  figure
  bw = cfg.bw;

  if bw
    % convert to greyscale image
    img = mean(img, 3);
    imagesc(img);
    colormap gray
  else
    % plot as RGB image
    image(img);
  end

  hold on
  axis equal
  axis off
  axis xy

  % get the electrode positions
  pos = zeros(0,2);
  electrodehelp = [ ...
    '-----------------------------------------------------\n' ...
    'specify electrode locations\n' ...
    'press the right mouse button to add another electrode\n' ...
    'press backspace on the keyboard to remove the last electrode\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(electrodehelp)
    disp(round(pos)); % values are integers/pixels
    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end
    
    switch k
      case 1
        pos = cat(1, pos, [x y]);
        % add it to the figure
        plot(x, y, 'b.');
        plot(x, y, 'yo');
      
      case 8
        if size(pos,1)>0
          % remove the last point
          pos = pos(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
        end

      case 'q'
        again = 0;
      
      otherwise
        warning('invalid button (%d)', k);
    end
  end

  % get the mask outline
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  maskhelp = [ ...
    '------------------------------------------------------------------------\n' ...
    'specify polygons for masking the topgraphic interpolation\n' ...
    'press the right mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(maskhelp);
    fprintf('\n');
    for i=1:length(polygon)
      fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
    end

    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end

    switch k
      case 1
        polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
        % add the last line segment to the figure
        if size(polygon{thispolygon},1)>1
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
        end
        plot(x, y, 'g.-');

      case 8 % backspace
        if size(polygon{thispolygon},1)>0
          % remove the last point
          polygon{thispolygon} = polygon{thispolygon}(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          % plot the electrode positions
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
          for i=1:length(polygon)
            x = polygon{i}(:,1);
            y = polygon{i}(:,2);
            if i~=thispolygon
              % close the polygon in the figure
              x(end) = x(1);
              y(end) = y(1);
            end
            plot(x, y, 'g.-');
          end
        end

      case 'c'
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end

      case 'q'
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
        end
        again = 0;

      otherwise
        warning('invalid button (%d)', k);
    end
  end
  % remember this set of polygons as the mask
  mask = polygon;


  % get the outline, e.g. head shape and sulci
  polygon = {};
  thispolygon = 1;
  polygon{thispolygon} = zeros(0,2);
  maskhelp = [ ...
    '-----------------------------------------------------------------------------------\n' ...
    'specify polygons for adding outlines (e.g. head shape and sulci) to the layout\n' ...
    'press the right mouse button to add another point to the current polygon\n' ...
    'press backspace on the keyboard to remove the last point\n' ...
    'press "c" on the keyboard to close this polygon and start with another\n' ...
    'press "n" on the keyboard to start with another without closing the current polygon\n' ...
    'press "q" on the keyboard to continue\n' ...
    ];
  again = 1;
  while again
    fprintf(maskhelp);
    fprintf('\n');
    for i=1:length(polygon)
      fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
    end

    try
      [x, y, k] = ginput(1);
    catch
      % this happens if the figure is closed
      return;
    end

    switch k
      case 1
        polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
        % add the last line segment to the figure
        if size(polygon{thispolygon},1)>1
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
        end
        plot(x, y, 'm.-');

      case 8 % backspace
        if size(polygon{thispolygon},1)>0
          % remove the last point
          polygon{thispolygon} = polygon{thispolygon}(1:end-1,:);
          % completely redraw the figure
          cla
          h = image(img);
          hold on
          axis equal
          axis off
          % plot the electrode positions
          plot(pos(:,1), pos(:,2), 'b.');
          plot(pos(:,1), pos(:,2), 'yo');
          for i=1:length(polygon)
            x = polygon{i}(:,1);
            y = polygon{i}(:,2);
            if i~=thispolygon
              % close the polygon in the figure
              x(end) = x(1);
              y(end) = y(1);
            end
            plot(x, y, 'm.-');
          end
        end

      case 'c'
        if size(polygon{thispolygon},1)>0
          x = polygon{thispolygon}(1,1);
          y = polygon{thispolygon}(1,2);
          polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
          % add the last line segment to the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'm.-');
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end

      case 'n'
        if size(polygon{thispolygon},1)>0
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
        end

      case 'q'
        again = 0;

      otherwise
        warning('invalid button (%d)', k);
    end
  end
  % remember this set of polygons as the outline
  outline = polygon;

  % convert electrode positions into a layout structure
  lay.pos = pos;
  nchans = size(pos,1);
  for i=1:nchans
    lay.label{i,1} = sprintf('chan%03d', i);
  end
  % add width and height for multiplotting
  d = dist(pos');
  for i=1:nchans
    d(i,i) = inf; % exclude the diagonal
  end
  mindist = min(d(:));
  lay.width  = ones(nchans,1) * mindist * 0.8;
  lay.height = ones(nchans,1) * mindist * 0.6;
  % add mask and outline polygons
  lay.mask    = mask;
  lay.outline = outline;

else
  fprintf('reverting to 151 channel CTF default\n');
  lay = readlay('CTF151.lay');
end

% FIXME there is a conflict between the use of cfg.style here and in topoplot
if ~strcmp(cfg.style, '3d')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check whether outline and mask are available
  % if not, add default "circle with triangle" to resemble the head
  % in case of "circle with triangle", the electrode positions should also be
  % scaled
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isfield(lay, 'outline') || ~isfield(lay, 'mask')
    rmax  = 0.5;
    l     = 0:2*pi/100:2*pi;
    HeadX = cos(l).*rmax;
    HeadY = sin(l).*rmax;
    NoseX = [0.18*rmax 0 -0.18*rmax];
    NoseY = [rmax-.004 rmax*1.15 rmax-.004];
    EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
    EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    % Scale the electrode positions to fit within a unit circle, i.e. electrode radius = 0.45
    ind_scale = strmatch('SCALE', lay.label);
    ind_comnt = strmatch('COMNT', lay.label);
    sel = setdiff(1:length(lay.label), [ind_scale ind_comnt]); % these are excluded for scaling
    x = lay.pos(sel,1);
    y = lay.pos(sel,2);
    xrange = range(x);
    yrange = range(y);
    % First scale the width and height of the box for multiplotting
    lay.width  = lay.width./xrange;
    lay.height = lay.height./yrange;
    % Then shift and scale the electrode positions
    lay.pos(:,1) = 0.9*((lay.pos(:,1)-min(x))/xrange-0.5);
    lay.pos(:,2) = 0.9*((lay.pos(:,2)-min(y))/yrange-0.5);
    % Define the outline of the head, ears and nose
    lay.outline{1} = [HeadX(:) HeadY(:)];
    lay.outline{2} = [NoseX(:) NoseY(:)];
    lay.outline{3} = [ EarX(:)  EarY(:)];
    lay.outline{4} = [-EarX(:)  EarY(:)];
    % Define the anatomical mask based on a circular head
    lay.mask{1} = [HeadX(:) HeadY(:)];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the montage, i.e. combine bipolar channels into a new representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(cfg.montage, 'no')
  Norg = length(cfg.montage.labelorg);
  Nnew = length(cfg.montage.labelnew);

  for i=1:Nnew
    cfg.montage.tra(i,:) = abs(cfg.montage.tra(i,:));
    cfg.montage.tra(i,:) = cfg.montage.tra(i,:) ./ sum(cfg.montage.tra(i,:));
  end
  % pretend it is a sensor structure, this achieves averaging after channel matching
  tmp.tra   = lay.pos;
  tmp.label = lay.label;
  new = apply_montage(tmp, cfg.montage);
  lay.pos   = new.tra;
  lay.label = new.label;
  % do the same for the width and height
  tmp.tra = lay.width(:);
  new = apply_montage(tmp, cfg.montage);
  lay.width = new.tra;
  tmp.tra = lay.height(:);
  new = apply_montage(tmp, cfg.montage);
  lay.height = new.tra;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add axes positions for comments and scale information if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~any(strcmp('COMNT', lay.label)) && strcmpi(cfg.style, '2d') && ~skipcomnt
  % add a placeholder for the comment in the upper left corner
  lay.label{end+1}  = 'COMNT';
  lay.width(end+1)  = mean(lay.width);
  lay.height(end+1) = mean(lay.height);
  X                 = min(lay.pos(:,1));
  Y                 = max(lay.pos(:,2));
  Y                 = min(lay.pos(:,2));
  lay.pos(end+1,:)  = [X Y];
elseif any(strcmp('COMNT', lay.label)) && skipcomnt
  % remove the scale entry
  sel = find(strcmp('COMNT', lay.label));
  lay.label(sel) = [];
  lay.pos(sel,:) = [];
  lay.width(sel) = [];
  lay.height(sel) = [];
end

if ~any(strcmp('SCALE', lay.label)) && strcmpi(cfg.style, '2d') && ~skipscale
  % add a placeholder for the scale in the upper right corner
  lay.label{end+1}  = 'SCALE';
  lay.width(end+1)  = mean(lay.width);
  lay.height(end+1) = mean(lay.height);
  X                 = max(lay.pos(:,1));
  Y                 = max(lay.pos(:,2));
  Y                 = min(lay.pos(:,2));
  lay.pos(end+1,:)  = [X Y];
elseif any(strcmp('SCALE', lay.label)) && skipscale
  % remove the scale entry
  sel = find(strcmp('SCALE', lay.label));
  lay.label(sel) = [];
  lay.pos(sel,:) = [];
  lay.width(sel) = [];
  lay.height(sel) = [];
end

% to plot the layout for debugging, you can use this code snippet
if strcmp(cfg.feedback, 'yes') && strcmpi(cfg.style, '2d')
  tmpcfg = [];
  tmpcfg.layout = lay;
  layoutplot(tmpcfg);
end

% to write the layout to a text file, you can use this code snippet
if ~isempty(cfg.output) && strcmpi(cfg.style, '2d')
  fprintf('writing layout to ''%s''\n', cfg.output);
  fid = fopen(cfg.output, 'wt');
  for i=1:numel(lay.label)
    fprintf(fid, '%d %f %f %f %f %s\n', i, lay.pos(i,1), lay.pos(i,2), lay.width(i), lay.height(i), lay.label{i});
  end
  fclose(fid);
elseif ~isempty(cfg.output) && strcmpi(cfg.style, '3d')
  % the layout file format does not support 3D positions, furthermore for
  % a 3D layout the width and height are currently set to NaN
  error('writing a 3D layout to an output file is not supported');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% read the layout information from the ascii file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lay = readlay(filename)
if ~exist(filename, 'file')
  error(sprintf('could not open layout file: %s', filename));
end
[chNum,X,Y,Width,Height,Lbl,Rem] = textread(filename,'%f %f %f %f %f %q %q');

if length(Rem)<length(Lbl)
  Rem{length(Lbl)} = [];
end

for i=1:length(Lbl)
  if ~isempty(Rem{i})
    % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
    Lbl{i} = [Lbl{i} ' ' Rem{i}];
  end
end
lay.pos    = [X Y];
lay.width  = Width;
lay.height = Height;
lay.label  = Lbl;
return % function readlay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% convert 3D electrode positions into 2D layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lay = sens2lay(sens, rz, method, style)

fprintf('creating layout for %s system\n', senstype(sens));
% apply rotation
if isempty(rz)
  switch senstype(sens)
    case {'ctf151', 'ctf275', 'bti148', 'bti248', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar'}
      rz = 90;
    case {'neuromag122', 'neuromag306'}
      rz = 0;
    case 'electrode'
      rz = 90;
    otherwise
      rz = 0;
  end
end
sens.pnt = warp_apply(rotate([0 0 rz]), sens.pnt, 'homogenous');

% use helper function for 3D layout
[pnt, label] = channelposition(sens);

if strcmpi(style, '3d')
  lay.pos   = pnt;
  lay.label = label;
else
  prj = elproj(pnt, method);
  d = dist(prj');
  d(find(eye(size(d)))) = inf;
  mindist = min(d(:));
  X = prj(:,1);
  Y = prj(:,2);
  Width  = ones(size(X)) * mindist * 0.8;
  Height = ones(size(X)) * mindist * 0.6;
  lay.pos    = [X Y];
  lay.width  = Width;
  lay.height = Height;
  lay.label  = label;
end
