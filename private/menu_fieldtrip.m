function menu_fieldtrip(handle, cfg, allowsubplot)

% MENU_FIELDTRIP adds a FieldTrip-specific menu to a figure.
%
% See also MENU_VIEWPOINT

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
    uimenu(ftmenu, 'Label', 'Show pipeline', 'Callback', {@cb_pipeline, cfg});
  end
  uimenu(ftmenu, 'Label', 'About', 'Callback', @cb_about);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_pipeline(handle, eventdata, varargin)
ft_analysispipeline([], varargin{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_about(handle, eventdata, varargin)
msgbox( {
  'FieldTrip is the MATLAB toolbox for MEG and EEG analysis that is being developed at the Centre for Cognitive Neuroimaging of the Donders Institute for Brain, Cognition and Behaviour together with collaborating institutes. The FieldTrip software is released as open source under the GNU general public license.'
  ''
  sprintf('This is FieldTrip, version %s', ft_version)
  } );
