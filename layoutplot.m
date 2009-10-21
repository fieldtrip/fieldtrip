function layoutplot(cfg, data)

% LAYOUTPLOT makes a figure with the 2-D layout of the channel positions
% for topoplotting and the individual channel axes (i.e. width and height
% of the subfigures) for multiplotting. A correct 2-D layout is a
% prerequisite  for plotting the topographical distribution of the
% potential or field distribution, or for plotting timecourses in a
% topographical arrangement.
%
% This function uses the same configuration options as prepare_layout and
% as the topoplotting and multiplotting functions. The difference is that
% this function plots the layout without any data, which facilitates
% the validation of your 2-D layout.
%
% Use as
%   layoutplot(cfg, data)
%
% There are several ways in which a 2-D layout can be made: it can be read
% directly from a *.lay file, it can be created based on 3-D electrode or
% gradiometer positions in the configuration or in the data, or it can be
% created based on the specification of an electrode of gradiometer file.
%
% You can specify either one of the following configuration options
%   cfg.layout      filename containg the layout
%   cfg.rotate      number, rotation around the z-axis in degrees (default = [], which means automatic)
%   cfg.projection  string, 2D projection method can be 'stereographic', 'ortographic', 'polar', 'gnomic' or 'inverse' (default = 'orthographic')
%   cfg.elec        structure with electrode positions, or
%   cfg.elecfile    filename containing electrode positions
%   cfg.grad        structure with gradiometer definition, or
%   cfg.gradfile    filename containing gradiometer definition
%   cfg.output      filename to which the layout will be written (default = [])
%   cfg.montage     'no' or a montage structure (default = 'no')
%   cfg.image       filename, use an image to construct a layout (e.g. usefull for ECoG grids)
%
% Alternatively the layout can be constructed from either
%   data.elec     structure with electrode positions
%   data.grad     structure with gradiometer definition
%
% Alternatively, you can specify
%   cfg.layout = 'ordered'
% which will give you a 2-D ordered layout. Note that this is only suited
% for multiplotting and not for topoplotting.
%
% See also prepare_layout, topoplotER, topoplotTFR, multiplotER, multiplotTFR

% Undocumented options
%   cfg.montage

% Copyright (C) 2006-2008, Robert Oostenveld
%
% $Log: layoutplot.m,v $
% Revision 1.11  2009/09/10 12:30:00  roboos
% added option for plotting an arrow for each of the bipolar electrode pairs, requires a montage and a layout
%
% Revision 1.10  2009/05/12 12:35:15  roboos
% moved plotting to seperate function (plot_lay) where it can be reused
%
% Revision 1.9  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.8  2008/09/22 12:44:46  roboos
% added support for cfg.image, some small changes to documentation and to plotting style
%
% Revision 1.7  2008/07/16 09:27:36  roboos
% added support for creating a layout based on a figure, including mask and outline for topoplot
%
% Revision 1.6  2007/03/21 14:47:06  chrhes
% updated some comments
%
% Revision 1.5  2007/03/20 10:44:37  roboos
% updated documentation, change figure axis, changed whitespaces
%
% Revision 1.4  2007/03/20 09:46:31  chrhes
% added some checks to determine whether the field cfg.layout already contains a
% valid layout structure that can be used, thereby avoiding a redundant call to
% the function PREPARE_LAYOUT
%
% Revision 1.3  2007/03/14 08:56:00  roboos
% changed some documentation
%
% Revision 1.2  2007/03/14 08:53:29  roboos
% changed from using private/createlayout to prepare_layout, changed the lay
% structure, added some help to clarify the usefullness of this function
%
% Revision 1.1  2006/06/01 11:51:45  roboos
% new implementation
%

fieldtripdefs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic check/initialization of input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1) || (nargin>2), error('incorrect number of input arguments'); end;
if (nargin<2), data = []; end;

if ~isstruct(cfg) && ~isempty(cfg), error('argument cfg must be a structure'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract/generate layout information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lay = [];

% try to use the layout structure from input if specified
if isfield(cfg, 'layout')
  % brief check to determine if cfg.layout is a valid layout (lay) structre
  if isstruct(cfg.layout)
    if all(isfield(cfg.layout, {'pos';'width';'height';'label'}))
      lay = cfg.layout;
    end
  end
end

% otherwise create the layout structure
if isempty(lay), lay = prepare_layout(cfg, data); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all details pertaining to the layout in one figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

if isfield(cfg, 'image') && ~isempty(cfg.image)
  % start with the background image
  fprintf('reading background image from %s\n', cfg.image);
  img = imread(cfg.image);
  img = flipdim(img, 1); % in combination with "axis xy"

  bw = 1;

  if bw
    % convert to greyscale image
    img = mean(img, 3);
    imagesc(img);
    colormap gray
  else
    % plot as RGB image
    image(img);
  end
  axis equal
  axis off
  axis xy
end

plot_lay(lay, 'point', true, 'box', true, 'label', true, 'mask', true, 'outline', true);

% the following code can be used to verify a bipolar montage, given the
% layout of the monopolar channels 
if isfield(cfg, 'montage') && ~isempty(cfg.montage)
  fprintf('plotting an arrow for each of the bipolar electrode pairs\n');
  % the arrow begins at the +1 electrode
  % the arrow ends   at the -1 electrode
  for i=1:length(cfg.montage.labelnew)
    begindx = find(cfg.montage.tra(i,:)==+1);
    endindx = find(cfg.montage.tra(i,:)==-1);
    if ~numel(begindx)==1 || ~numel(endindx)==1
      % the re-referenced channel does not seem to be a bipolar pair
      continue
    end
    % find the position of the begin and end of the arrow
    beglab = cfg.montage.labelorg{begindx};
    endlab = cfg.montage.labelorg{endindx};
    begindx = find(strcmp(lay.label, beglab)); % the index in the layout
    endindx = find(strcmp(lay.label, endlab)); % the index in the layout
    if ~numel(begindx)==1 || ~numel(endindx)==1
      % one of the channels in the bipolar pair does not seem to be in the layout
      continue
    end

    begpos = lay.pos(begindx,:);
    endpos = lay.pos(endindx,:);
    arrow(begpos, endpos, 'Length', 5)

  end % for all re-referenced channels
end % if montage


