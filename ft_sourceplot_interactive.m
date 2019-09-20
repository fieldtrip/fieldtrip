function ft_sourceplot_interactive(cfg, varargin)

% FT_SOURCEPLOT_INTERACTIVE provides a rapid way to plot 3D surface
% renderings of pos_time or pos_freq functional data, and interactively
% explore them. One figure is created with surface plots of the individual
% conditions, and by default a plot of the functional data averaged over
% the entire cortex is created over time (or frequency). Users can click in
% the line graph to shift the time point for which the functional data is
% shown in the surface plots. Additionally, users can Shift+Click in the
% surface plots to add a "virtual electrode", for which a new line graph
% figure will be created.
%
% Input data needs to be source+mesh, so has to contain a tri, pos, and one
% functional field plus a time- or frequency axis.
%
% Use e.g. like so:
%
% cfg = [];
% cfg.data_labels = {'Congruent', 'Incongruent'};
% ft_sourceplot_interactive(cfg, sourceFC, sourceFIC);
%
% Copyright (C) 2019 Eelke Spaak, Donders Institute. e.spaak@donders.ru.nl

% these are used by the ft_preamble/ft_postamble function and scripts
ft_nargin = nargin;
ft_nargout = nargout;

% call the set of 'standard' preambles (see ft_examplefunction for details)
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return;
end

% validate the input
if numel(varargin) < 1
  ft_error('this function requires at least one data input argument');
end
for k = 1:numel(varargin)
  varargin{k} = ft_checkdata(varargin{k}, 'datatype', {'source+mesh'}, 'feedback', 'yes');
  if k > 1 && (~isequaln(varargin{k}.pos, varargin{1}.pos) || ...
    ~isequaln(varargin{k}.tri, varargin{1}.tri))
    % ensure all input arguments are expressed on the same mesh
    ft_error('input arguments for plotting need to all have identical .pos and .tri');
  end
end

if isfield(varargin{1}, 'time') || isfield(varargin{1}, 'freq')
  % make a selection of the time and/or frequency dimension
  tmpcfg = keepfields(cfg, {'frequency', 'avgoverfreq', 'keepfreqdim', 'latency',...
    'avgovertime', 'keeptimedim', 'showcallinfo'});
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
end

% validate the cfg options
cfg.parameter = ft_getopt(cfg, 'parameter', 'pow');

% check whether we're dealing with time or frequency data
dimord = getdimord(varargin{1}, cfg.parameter);
if ~ismember(dimord, {'pos_time', 'pos_freq'})
  ft_error('functional data must be pos_time or pos_freq');
end

% set a sensible x-axis label for the plots over time or freq
cfg.time_label = ft_getopt(cfg, 'time_label', []);
if isempty(cfg.time_label)
  if strcmp(dimord, 'pos_time')
    cfg.time_label = 'Time (s)';
  elseif strcmp(dimord, 'pos_freq')
    cfg.time_label = 'Frequency (Hz)';
  end
end

if strcmp(dimord, 'pos_time')
  xdat = varargin{1}.time;
elseif strcmp(dimord, 'pos_freq')
  xdat = varargin{1}.freq;
end

% optionally load an atlas
if isfield(cfg, 'atlas') && ischar(cfg.atlas)
  cfg.atlas = ft_read_atlas(cfg.atlas);
end

% other defaults are set in the lower-level object

% fetch the functional data
data = cellfun(@(x) x.(cfg.parameter), varargin, 'uniformoutput', false);

% set up the arguments
keyval = struct2keyval(cfg);
keyval = [keyval {'tri', varargin{1}.tri, 'pos', varargin{1}.pos, 'data', data, 'time', xdat}];

% and launch the viewer
viewer = ft_interactivesourceviewer(keyval{:});
viewer.show();

ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble savefig

end
