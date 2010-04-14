function ft_layoutplot(cfg, data)

% FT_LAYOUTPLOT makes a figure with the 2-D layout of the channel positions
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
%   ft_layoutplot(cfg, data)
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
% See also ft_prepare_layout, ft_topoplotER, ft_topoplotTFR, ft_multiplotER, ft_multiplotTFR

% Undocumented options
%   cfg.montage

% Copyright (C) 2006-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
if isempty(lay), lay = ft_prepare_layout(cfg, data); end;

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


