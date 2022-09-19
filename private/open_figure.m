function h = open_figure(cfg)

% OPEN_FIGURE is a helper function to open a figure with some specific settings
% consistent over all FieldTrip functions that do plotting and/or that show a
% graphical user interface.
%
% See also GCA, GCF, GROOT, 

cfg.figure      = ft_getopt(cfg, 'figure', 'yes');
cfg.visible     = ft_getopt(cfg, 'visible', 'yes');
cfg.position    = ft_getopt(cfg, 'position'); % use the default position
cfg.renderer    = ft_getopt(cfg, 'renderer'); % let MATLAB decide on the default
cfg.figurename  = ft_getopt(cfg, 'figurename');
cfg.title       = ft_getopt(cfg, 'title');

switch cfg.figure
  case {'new', 'yes'}
    figopt = {};
    if ~istrue(cfg.visible)
      figopt = ft_setopt(figopt, 'visible', 'off');
    end
    if ~isempty(cfg.position)
      figopt = ft_setopt(figopt, 'position', cfg.position);
    end
    
    % check whether there is already a figure open
    h = get(0, 'CurrentFigure');

    if isempty(h)
      % there is no figure open yet, make a new one
      h = figure(figopt{:});
    elseif isempty(get(h, 'Children'))
      % there is an empty figure open, use that and update it
      if ~isempty(figopt)
        set(h, figopt{:});
      end
    else
      % there is another figure open, make a new one
      h = figure(figopt{:});
    end
    
  case {'clf', 'no'}
    % use the current figure and clear it
    % this will open a new one if there is no figure yet
    h = clf;

  case 'gca'
    % use the current axes but do not clear them
    % this will open a new figure if there is no figure yet
    h = gca;

  case 'gcf'
    % use the current figure but do not clear it
    % this will open a new figure if there is no figure yet
    h = gcf;
  
  otherwise
    % assume that it specifies a figure handle
    h = cfg.figure;
end

assert(ishandle(h), 'failed to open figure');

if ~isempty(cfg.figurename)
  % this appears as the name of the window
  set(h, 'name', cfg.figurename);
end

if ~isempty(cfg.renderer)
  % set figure renderer
  set(h, 'Renderer', cfg.renderer);
end

if ~isempty(cfg.title)
  % this appears above the axes
  title(cfg.title);
end
