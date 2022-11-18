function [varargout] = bivariate_common(cfg, varargin)

% BIVARIATE_COMMON makes a selection for a specific reference channel from a
% bivariate (i.e. connectivity) dataset and returns that selection as a univariate
% dataset. This is used in singleplot/multiplot/topoplot for both ER and TFR data.
%
% Use as
%   [varargout] = bivariate_common(cfg, varargin)
%
% See also TOPOPLOT_COMMON

% reference channel is required
if ~isfield(cfg, 'refchannel') || isempty(cfg.refchannel)
  ft_error('no reference channel is specified');
end

% check for refchannel being part of selection
if ~strcmp(cfg.refchannel, 'gui')
  if isfield(varargin{1}, 'labelcmb')
    cfg.refchannel = ft_channelselection(cfg.refchannel, unique(varargin{1}.labelcmb(:)));
  else
    cfg.refchannel = ft_channelselection(cfg.refchannel, varargin{1}.label);
  end
  if (isfield(varargin{1}, 'label')    && ~any(ismember(varargin{1}.label,       cfg.refchannel))) || ...
      (isfield(varargin{1}, 'labelcmb') && ~any(ismember(varargin{1}.labelcmb(:), cfg.refchannel)))
    ft_error('cfg.refchannel is a not present in the (selected) channels)')
  end
end

% Interactively select the reference channel
if strcmp(cfg.refchannel, 'gui')
  tmpcfg = keepfields(cfg, {'channel', 'layout', 'commentpos', 'scalepos', 'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  lay = ft_prepare_layout(tmpcfg, varargin{1});
  % Open a single figure with the channel layout, the user can click on a reference channel
  h = clf;
  ft_plot_layout(lay, 'box', false);
  title('Select the reference channel by dragging a selection window, more than 1 channel can be selected...');
  % add the channel information to the figure
  info       = guidata(gcf);
  info.x     = lay.pos(:, 1);
  info.y     = lay.pos(:, 2);
  info.label = lay.label;
  info.dataname = '';
  guidata(h, info);
  set(h, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@make_selection, cfg, varargin{:}}, 'event', 'WindowButtonUpFcn'});
  set(h, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@make_selection, cfg, varargin{:}}, 'event', 'WindowButtonDownFcn'});
  set(h, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@make_selection, cfg, varargin{:}}, 'event', 'WindowButtonMotionFcn'});
else
  
  make_figure(cfg, varargin{~cellfun(@isnumeric,varargin)});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which creates the figure with the selected reference channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_figure(cfg, varargin)

if isfield(cfg, 'inputfile')
  % the reading has already been done and varargin contains the data
  cfg = rmfield(cfg, 'inputfile');
end

for i=1:numel(varargin)
  dimord = getdimord(varargin{i}, cfg.parameter);
  dimtok = tokenize(dimord, '_');
  
  if isequal(dimtok{1}, 'chancmb')
    % convert 2-dimensional channel matrix to a single dimension
    if isempty(cfg.directionality)
      sel1 = match_str(varargin{i}.labelcmb(:, 1), cfg.refchannel);
      sel2 = match_str(varargin{i}.labelcmb(:, 2), cfg.refchannel);
    elseif strcmp(cfg.directionality, 'outflow')
      sel1 = match_str(varargin{i}.labelcmb(:, 1), cfg.refchannel);
      sel2 = [];
    elseif strcmp(cfg.directionality, 'inflow')
      sel1 = [];
      sel2 = match_str(varargin{i}.labelcmb(:, 2), cfg.refchannel);
    end
    fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.parameter);
    if length(sel1)+length(sel2)==0
      ft_error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
    end
    
    % combine the selection in the first and second column
    selindx = [sel1; sel2];
    sellab  = [varargin{i}.labelcmb(sel1, 2); varargin{i}.labelcmb(sel2, 1)];
    [sellab, tmp] = unique(sellab, 'stable');
    selindx = selindx(tmp);
    
    % take the selected rows from the data and the corresponding labels
    varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)(selindx, :, :);
    varargin{i}.label           = sellab;
    varargin{i} = rmfield(varargin{i}, 'labelcmb');
    
    % update the dimord
    dimtok{1} = 'chan';
    varargin{i} = removefields(varargin{i}, {'dimord', [cfg.parameter 'dimord']});
    varargin{i}.dimord = sprintf('%s_', dimtok{1:end});
    varargin{i}.dimord(end) = []; % remove the trailing "_"
    
  elseif isequal(dimtok{1}, 'chan') && isequal(dimtok{2}, 'chan')
    sel = match_str(varargin{i}.label, cfg.refchannel);
    siz = [size(varargin{i}.(cfg.parameter)) 1];
    
    if strcmp(cfg.directionality, 'inflow') || isempty(cfg.directionality)
      % The interpretation of 'inflow' and 'outflow' depend on the definition in the
      % bivariate representation of the data. In FieldTrip the row index 'causes' the
      % column index channel
      sel1    = 1:siz(1);
      sel2    = sel;
      meandir = 2;
      
    elseif strcmp(cfg.directionality, 'outflow')
      sel1    = sel;
      sel2    = 1:siz(1);
      meandir = 1;
      
    elseif strcmp(cfg.directionality, 'ff-fd')
      ft_error('cfg.directionality = ''ff-fd'' is not supported anymore, you have to manually subtract the two prior to plotting');
      
    elseif strcmp(cfg.directionality, 'fd-ff')
      ft_error('cfg.directionality = ''fd-ff'' is not supported anymore, you have to manually subtract the two prior to plotting');
    end % if directionality
    
    % Make a univariate selection of the data
    varargin{i}.(cfg.parameter) = mean(varargin{i}.(cfg.parameter)(sel1, sel2,:,:), meandir);
    siz(meandir) = [];
    varargin{i}.(cfg.parameter) = reshape(varargin{i}.(cfg.parameter), siz);
    
    % Update the dimord
    varargin{i} = removefields(varargin{i}, {'dimord', [cfg.parameter 'dimord']});
    varargin{i}.dimord = sprintf('%s_', dimtok{2:end});
    varargin{i}.dimord(end) = []; % remove the trailing "_"
    
  else
    error('unexpected dimord');
    
  end % if sparse or full
end % for varargin

% Remove these fields from the configuration
fn = {'originalfunction', 'inputfile', 'refchannel'};

% This applies to the topoplots
cfg.highlight = 'on';
cfg.highlightsymbol  = '.';
cfg.highlightcolor   = 'r';
cfg.highlightsize = 20;
cfg.highlightchannel =  cfg.refchannel;

switch cfg.originalfunction
  case 'ft_multiplotER'
    ft_multiplotER(removefields(cfg, fn), varargin{:});
  case 'ft_multiplotTFR'
    ft_multiplotTFR(removefields(cfg, fn), varargin{:});
  case 'ft_singleplotER'
    ft_singleplotER(removefields(cfg, fn), varargin{:});
  case 'ft_singleplotTFR'
    ft_singleplotTFR(removefields(cfg, fn), varargin{:});
  case 'ft_topoplotER'
    ft_topoplotER(removefields(cfg, fn), varargin{:});
  case 'ft_topoplotTFR'
    ft_topoplotTFR(removefields(cfg, fn), varargin{:});
  otherwise
    error('unsupported cfg.originalfunction')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called by ft_select_channel in case cfg.refchannel='gui'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_selection(label, cfg, varargin)
fprintf('selected cfg.refchannel = ''%s''\n', join_str(', ', cfg.refchannel));
cfg.refchannel = label;
make_figure(cfg, varargin{:})
