function [cfg] = ft_connectivityplot(cfg, varargin)

% FT_CONNECTIVITYPLOT plots channel-level frequency resolved connectivity. The
% data are rendered in a square grid of subplots, each subplot containing the
% connectivity spectrum between the two respective channels.
%
% Use as
%   ft_connectivityplot(cfg, data)
% where the first input argument is a configuration structure (see below)
% and the input data is a structure obtained from  FT_CONNECTIVITYANALYSIS
% using a frequency-domain connectivity metric. Consequently the input data
% should have a dimord of 'chan_chan_freq', or 'chan_chan_freq_time'.
%
% The configuration can have the following options
%   cfg.parameter   = string, the functional parameter to be plotted (default = 'cohspctrm')
%   cfg.xlim        = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [xmin xmax] (default = 'maxmin')
%   cfg.ylim        = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.zlim        = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.channel     = list of channels to be included for the plotting (default = 'all'), see FT_CHANNELSELECTION for details
%
% See also FT_CONNECTIVITYANALYSIS, FT_CONNECTIVITYSIMULATION, FT_MULTIPLOTCC, FT_TOPOPLOTCC

% Copyright (C) 2011-2017, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'zparam',     'parameter'});
cfg = ft_checkconfig(cfg, 'renamed',    {'color',      'linecolor'});
cfg = ft_checkconfig(cfg, 'renamed',    {'graphcolor', 'linecolor'});

% set the defaults
cfg.channel     = ft_getopt(cfg, 'channel',   'all');
cfg.parameter   = ft_getopt(cfg, 'parameter', 'cohspctrm');
cfg.zlim        = ft_getopt(cfg, 'zlim',      'maxmin');
cfg.ylim        = ft_getopt(cfg, 'ylim',      'maxmin');
cfg.xlim        = ft_getopt(cfg, 'xlim',      'maxmin');
cfg.linecolor   = ft_getopt(cfg, 'linecolor', 'brgkywrgbkywrgbkywrgbkyw');
cfg.linestyle   = ft_getopt(cfg, 'linestyle', '-');
cfg.linewidth   = ft_getopt(cfg, 'linewidth', 0.5);
cfg.visible     = ft_getopt(cfg, 'visible',   'on');
cfg.renderer    = ft_getopt(cfg, 'renderer',  []); % let MATLAB decide on the default

% check if the input data is valid for this function
% ensure that the input is correct
Ndata = numel(varargin);
dtype = cell(Ndata, 1);
iname = cell(Ndata+1, 1);
for k = 1:Ndata
  if ischar(cfg.parameter)
    cfg.parameter = repmat({cfg.parameter}, [1 Ndata]);
  end

  % check whether all requested parameters are the same. If not, rename
  % this, because otherwise a call to ft_selectdata (below) won't work
  if ~all(strcmp(cfg.parameter,cfg.parameter{1}))
    fprintf('different types of connectivity are to be displayed in the same figure\n');
    varargin{k}.connectivity = varargin{k}.(cfg.parameter{k});
    varargin{k} = rmfield(varargin{k}, cfg.parameter{k});
    if isfield(varargin{k}, [cfg.parameter{k} 'dimord'])
      varargin{k}.connectivitydimord = varargin{k}.([cfg.parameter{k} 'dimord']);
      varargin{k} = rmfield(varargin{k}, [cfg.parameter{k} 'dimord']);
    end
    cfg.parameter{k} = 'connectivity';
  else
    % don't worry
  end

  % check if the input data is valid for this function
  varargin{k} = ft_checkdata(varargin{k}, 'datatype', {'timelock', 'freq'});
  dtype{k}    = ft_datatype(varargin{k});

  % convert into the the supported dimord
  switch varargin{k}.dimord
    case {'chan_chan_freq' 'chan_chan_freq_time'}
      % that's ok
    case {'chancmb_freq' 'chancmb_freq_time'}
      % convert into 'chan_chan_freq'
      varargin{k} = ft_checkdata(varargin{k}, 'cmbstyle', 'full');
    otherwise
      ft_error('the data should have a dimord of %s or %s', 'chan_chan_freq', 'chancmb_freq');
  end

  % this is needed for correct treatment of linecolor later on
  if nargin>1
    if ~isempty(inputname(k+1))
      iname{k+1} = inputname(k+1);
    else
      iname{k+1} = ['input',num2str(k, '%02d')];
    end
  else
    % not yet supported
    iname{k+1} = cfg.inputfile{k};
  end
end

if Ndata >1
  if ~all(strcmp(dtype{1}, dtype))
    ft_error('input data are of different type; this is not supported');
  end
end

if ischar(cfg.linecolor)
  % ensure it is a column vector of the right length
  cfg.linecolor = repmat(cfg.linecolor(:), ceil(Ndata/length(cfg.linecolor)), 1);
  cfg.linecolor = cfg.linecolor(1:Ndata);
elseif isnumeric(cfg.linecolor)
  % ensure it is a Nx3 matrix of the right length
  cfg.linecolor = repmat(cfg.linecolor, ceil(Ndata/length(cfg.linecolor)), 1);
  cfg.linecolor = cfg.linecolor(1:Ndata,:);
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'channel', 'showcallinfo', 'trackcallinfo', 'trackconfig', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

% check presence of time / freq axes
hasfreq = isfield(varargin{1}, 'freq');
hastime = isfield(varargin{1}, 'time');
if hasfreq && hastime
  if Ndata>1
    ft_error('when the input data contains time-frequency representations, only a single data argument is allowed');
  end
  xparam = 'time';
  yparam = 'freq';
elseif hasfreq
  xparam = 'freq';
  yparam = '';
elseif hastime
  xparam = 'time';
  yparam = '';
end

% Get physical min/max range of x, i.e. time
if ~isnumeric(cfg.xlim)
  % Find maxmin throughout all varargins
  xmin = +inf;
  xmax = -inf;
  for i=1:Ndata
    xmin = nanmin([xmin varargin{i}.(xparam)]);
    xmax = nanmax([xmax varargin{i}.(xparam)]);
  end
  switch cfg.xlim
    case 'maxmin'
      % keep them as they are
    case 'maxabs'
      xmax = max(abs(xmax), abs(xmin));
      xmin = -xmax;
    case 'zeromax'
      xmin = 0;
    case 'minzero'
      xmax = 0;
    otherwise
      ft_error('invalid specification of cfg.xlim');
  end % switch
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end
cfg.xlim = [xmin xmax];

% Get physical min/max range of y, i.e. freq
if isempty(yparam)
  ymin = [];
  ymax = [];
elseif ~isnumeric(cfg.ylim)
  % Find maxmin throughout all varargins
  ymin = +inf;
  ymax = -inf;
  for i=1:Ndata
    ymin = nanmin([ymin varargin{i}.(yparam)]);
    ymax = nanmax([ymax varargin{i}.(yparam)]);
  end
  switch cfg.ylim
    case 'maxmin'
      % keep them as they are
    case 'maxabs'
      ymax = max(abs(ymax), abs(ymin));
      ymin = -ymax;
    case 'zeromax'
      ymin = 0;
    case 'minzero'
      ymax = 0;
    otherwise
      ft_error('invalid specification of cfg.ylim');
  end % switch
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end  
cfg.ylim = [ymin ymax];

% Get physical min/max range of z, which is the functional data:
if ~isnumeric(cfg.zlim)
  zmin = +inf;
  zmax = -inf;
  for k = 1:Ndata
    zmin = nanmin([zmin min(varargin{k}.(cfg.parameter{k})(:))]);
    zmax = nanmax([zmax max(varargin{k}.(cfg.parameter{k})(:))]);
  end
  switch cfg.zlim
    case 'maxmin'
      % keep them as they are
    case 'maxabs'
      zmax = max(abs(zmax), abs(zmin));
      zmin = -zmax;
    case 'zeromax'
      zmin = 0;
    case 'minzero'
      zmax = 0;
    otherwise
      ft_error('invalid specification of cfg.ylim');
  end % switch
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end
cfg.zlim = [zmin zmax];

% make the function recursive if Ndata>1
if Ndata>1
  data = varargin{1};
  tmpcfg = cfg;
  if ischar(cfg.parameter)
    % do nothing
  elseif iscell(cfg.parameter)
    tmpcfg.parameter = cfg.parameter{1};
  end
  ft_connectivityplot(tmpcfg, data);
  tmpcfg = cfg;

  if ischar(cfg.linecolor),        colorLabels = [iname{2} '='         tmpcfg.linecolor(1)     '\n'];
  elseif isnumeric(cfg.linecolor), colorLabels = [iname{2} '=' num2str(tmpcfg.linecolor(1, :)) '\n'];
  end

  for k = 2:Ndata
    if ischar(cfg.linecolor),     tmpcfg.linecolor = tmpcfg.linecolor(2:end);
    else isnumeric(cfg.linecolor),tmpcfg.linecolor = tmpcfg.linecolor(2:end,:);
    end

    tmpcfg.holdfig = 1;
    if ischar(cfg.parameter)
      % do nothing
    elseif iscell(cfg.parameter)
      tmpcfg.parameter = cfg.parameter{k};
    end
    ft_connectivityplot(tmpcfg, varargin{k});

    if ischar(cfg.linecolor);        colorLabels = [colorLabels iname{k+1} '='         tmpcfg.linecolor(1)     '\n'];
    elseif isnumeric(cfg.linecolor); colorLabels = [colorLabels iname{k+1} '=' num2str(tmpcfg.linecolor(1, :)) '\n'];
    end
  end

  ft_plot_text(0.5, (numel(varargin{k}.label)+1).*1.2-0.5, sprintf(colorLabels), 'horizontalalignment', 'right');

  return;
else
  data = varargin{1};
end

if ~isfield(data, cfg.parameter{1})
  ft_error('the data does not contain the requested parameter %s', cfg.parameter{1});
end

% get the selection of the data
tmpcfg           = [];
if hasfreq && hastime
  tmpcfg.latency   = cfg.xlim;
  tmpcfg.frequency = cfg.ylim;
elseif hasfreq
  tmpcfg.frequency = cfg.xlim;
elseif hastime
  tmpcfg.latency   = cfg.xlim;
end
data             = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

dat   = data.(cfg.parameter{k});
nchan = numel(data.label);
if hasfreq, nfreq = numel(data.freq); end
if hastime, ntime = numel(data.time); end

if (isfield(cfg, 'holdfig') && cfg.holdfig==0) || ~isfield(cfg, 'holdfig')
  cla;
  hold on;
end

for k = 1:nchan
  for m = 1:nchan
    if k~=m
      ix  = k;
      iy  = nchan - m + 1;
      % use the convention of the row-channel causing the column-channel
      if hastime && hasfreq && ntime>1
        tmp = reshape(dat(m,k,:,:), [nfreq ntime]);
        ft_plot_matrix(tmp, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'clim', cfg.zlim, 'box', 'yes');
      elseif hasfreq
        tmp = reshape(dat(m,k,:), [nfreq 1]);
        ft_plot_vector(tmp, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'vlim', cfg.zlim, 'box', 'yes', 'color', cfg.linecolor(1,:), 'linewidth', cfg.linewidth, 'style', cfg.linestyle);
      elseif hastime
        ft_error('plotting data with only a time axis is not supported yet');
      end
      if k==1
        % first column, plot scale on y axis
        if hastime && hasfreq && ntime>1
          val1 = cfg.ylim(1);
          val2 = cfg.ylim(2);
        elseif hasfreq
          val1 = cfg.zlim(1);
          val2 = cfg.zlim(2);
        elseif hastime
        end
        fontsize = 10;
        ft_plot_text(ix.*1.2-0.5, iy.*1.2-0.5, num2str(val1,3), 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle', 'FontSize', fontsize);
        ft_plot_text(ix.*1.2-0.5, iy.*1.2+0.5, num2str(val2,3), 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle', 'FontSize', fontsize);
      end
      if m==nchan
        % bottom row, plot scale on x axis
        fontsize = 10;
        ft_plot_text(ix.*1.2-0.5, iy.*1.2-0.5, num2str(cfg.xlim(1),3), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'FontSize', fontsize);
        ft_plot_text(ix.*1.2+0.5, iy.*1.2-0.5, num2str(cfg.xlim(2),3), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'FontSize', fontsize);
      end
    end
  end
end

% add channel labels on grand X and Y axes
for k = 1:nchan
  ft_plot_text(0.5,     (nchan + 1 - k).*1.2,     data.label{k}, 'horizontalalignment', 'right');
  ft_plot_text(k.*1.2,  (nchan + 1)    .*1.2-0.5, data.label{k}, 'horizontalalignment', 'left', 'rotation', 90);
end

% add 'from' and 'to' labels
ft_plot_text(-0.5,            (nchan + 1)/1.7,     '\it{from}', 'interpreter', 'tex', 'rotation', 90);
ft_plot_text((nchan + 1)/1.7, (nchan + 1)*1.2+0.4, '\it{to}', 'interpreter', 'tex');

axis([-0.2 (nchan+1).*1.2+0.2 0 (nchan+1).*1.2+0.2]);
axis off;

set(gcf, 'color', [1 1 1]);

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

% set the figure window title
if ~isempty(dataname)
  set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), mfilename, join_str(', ', dataname)));
else
  set(gcf, 'Name', sprintf('%d: %s', double(gcf), mfilename));
end
set(gcf, 'NumberTitle', 'off');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end
