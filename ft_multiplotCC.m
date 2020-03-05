function [cfg] = ft_multiplotCC(cfg, data)

% FT_MULTIPLOTCC visualises the coherence between channels by using
% multiple topoplots. The topoplot at a given channel location shows the
% coherence of that channel with all other channels.
%
% Use as
%   ft_multiplotCC(cfg, data)
%
% See also FT_PREPARE_LAYOUT, FT_TOPOPLOTCC, FT_CONNECTIVITYPLOT

% Undocumented local options:
% cfg.layout  = layout filename or a structure produced by prepare_layout
% cfg.xlim
% cfg.parameter
% This function requires input from FT_FREQSTATISTICS_SHIFTPREDICT
% This function should be rewritten, using the clean topoplot implementation

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen, Robert Oostenveld
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data);

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'xparam'});

% set the defaults
cfg.xlim      = ft_getopt(cfg, 'xlim', 'all');
cfg.parameter = ft_getopt(cfg, 'parameter', 'avg.icohspctrm');
cfg.renderer  = ft_getopt(cfg, 'renderer'); % let MATLAB decide on the default

if strcmp(cfg.parameter, 'avg.icohspctrm') && ~issubfield(data, 'avg.icohspctrm')
  data.avg.icohspctrm = abs(imag(data.avg.cohspctrm));
end

if strcmp(data.dimord, 'chan_chan_freq')
  % reshape input-data, such that ft_topoplotTFR will take it
  cnt = 1;
  siz = size(data.prob);
  data.labelcmb = cell(siz(1)*siz(2),2);
  data.prob = reshape(data.prob, [siz(1)*siz(2) siz(3)]);
  data.stat = reshape(data.stat, [siz(1)*siz(2) siz(3)]);
  for j = 1:length(data.label)
    for k = 1:length(data.label)
      data.labelcmb(cnt,:) = [data.label(k) data.label(j)];
      cnt = cnt + 1;
    end
  end
  tmpdata = data;
else
  dat   = getsubfield(data, cfg.parameter);
  scale = [0 max(dat(:))-0.2];
end

if isfield(cfg, 'xparam')
  xparam = getsubfield(data, cfg.xparam);
  if ~strcmp(cfg.xlim, 'all')
    fbin = [nearest(xparam, cfg.xlim(1)) nearest(xparam, cfg.xlim(2))];
  else
    fbin = [xparam(1) xparam(end)];
  end
end

% read or create the layout that will be used for plotting
tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'scalepos', 'elec', 'grad', 'opto', 'showcallinfo'});
layout = ft_prepare_layout(tmpcfg);
ft_plot_layout(layout, 'box', false, 'label', 'no', 'point', 'no');

X = layout.pos(:,1);
Y = layout.pos(:,2);
Width = layout.width;
Height = layout.height;
Lbl = layout.label;
chNum = numel(layout.label);

xScaleFac = 1/(max(Width)+ max(X) - min(X));
yScaleFac = 1/(max(Height)+ max(Y) - min(Y));

Xpos = xScaleFac*(X-min(X));
Ypos = 0.9*yScaleFac*(Y-min(Y));

for k=1:length(chNum) - 2
  subplotOL('position',[Xpos(k) Ypos(k)+(Height(k)*yScaleFac) Width(k)*xScaleFac*2 Height(k)*yScaleFac*2])
  tmpcfg = [];
  tmpcfg.layout = layout;
  if exist('tmpdata', 'var')
    tmpcfg.style      = 'straight';
    tmpcfg.marker     = 'off';
    try
      tmpcfg.refmarker = strmatch(Lbl(k), data.reflabel);
    catch
      tmpcfg.refmarker = strmatch(Lbl(k), data.label);
    end
    tmpcfg.interplimits = 'electrodes';
    if isfield(cfg, 'xparam')
      tmpcfg.xparam = cfg.xparam;
      tmpcfg.xlim   = xparam;
    else
      tmpcfg.xparam = 'freq';
      tmpcfg.xlim   = [k-0.5 k+0.5];
    end
    tmpcfg.parameter  = cfg.parameter;
    tmpcfg.refchannel = Lbl(k);
    tmpcfg.colorbar   = 'no';
    tmpcfg.zlim       = scale;
    tmpcfg.grid_scale = 30;
    ft_topoplotTFR(tmpcfg, data);
    drawnow
  end
end

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

% set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end

