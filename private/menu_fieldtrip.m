function menu_fieldtrip(handle, cfg, allowsubplot)

% MENU_FIELDTRIP adds a FieldTrip-specific menu to a figure.
%
% See also MENU_ABOUT, MENU_PIPELINE

if nargin<2
  cfg = [];
end

if nargin<3
  allowsubplot = true;
end

if ~allowsubplot && numel(findobj(gcf, 'type', 'axes')) > 1
  % do not add a menu if there are multiple axes, i.e. subplots
  % these may contain different data structures
  return
end

if ft_platform_supports('uimenu')
  % delete any possibly existing previous menu, this is safe because delete([]) does nothing
  delete(findobj(handle, 'type', 'uimenu', 'label', 'FieldTrip'));

  % add the menu with callbacks
  ftmenu = uimenu(handle, 'Label', 'FieldTrip');
  if ~isempty(cfg)
    uimenu(ftmenu, 'Label', 'Show pipeline', 'Callback', {@menu_pipeline, cfg});
  end
  uimenu(ftmenu, 'Label', 'About', 'Callback', @menu_about);
end
