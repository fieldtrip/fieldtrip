function [cfg] = ft_singleplotER(cfg, varargin)

% FT_SINGLEPLOTER plots the event-related fields or potentials of a single
% channel or the average over multiple channels. Multiple datasets can be
% overlayed.
%
% Use as
%   ft_singleplotER(cfg, data)
% or
%   ft_singleplotER(cfg, data1, data2, ..., datan)
%
% The data can be an erp/erf produced by FT_TIMELOCKANALYSIS, a power
% spectrum produced by FT_FREQANALYSIS or connectivity spectrum produced by
% FT_CONNECTIVITYANALYSIS.
%
% The configuration can have the following parameters:
%   cfg.parameter     = field to be plotted on y-axis (default depends on data.dimord)
%                       'avg', 'powspctrm' or 'cohspctrm'
%   cfg.maskparameter = field in the first dataset to be used for masking of data
%                       (not possible for mean over multiple channels, or when input contains multiple subjects
%                       or trials)
%   cfg.maskstyle     = style used for masking of data, 'box', 'thickness' or 'saturation' (default = 'box')
%   cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.ylim          = 'maxmin', 'maxabs', 'zeromax', 'minzero', or [ymin ymax] (default = 'maxmin')
%   cfg.channel       = nx1 cell-array with selection of channels (default = 'all'),
%                       see ft_channelselection for details
%   cfg.refchannel    = name of reference channel for visualising connectivity, can be 'gui'
%   cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see ft_timelockbaseline
%   cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
%   cfg.trials        = 'all' or a selection given as a 1xn vector (default = 'all')
%   cfg.fontsize      = font size of title (default = 8)
%   cfg.hotkeys       = enables hotkeys (up/down/left/right arrows) for dynamic x/y axis translation (Ctrl+) and zoom adjustment
%   cfg.interactive   = interactive plot 'yes' or 'no' (default = 'yes')
%                       in a interactive plot you can select areas and produce a new
%                       interactive plot when a selected area is clicked. multiple areas
%                       can be selected by holding down the shift key.
%   cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
%   cfg.linestyle     = linestyle/marker type, see options of the PLOT function (default = '-')
%                       can be a single style for all datasets, or a cell-array containing one style for each dataset
%   cfg.linewidth     = linewidth in points (default = 0.5)
%   cfg.graphcolor    = color(s) used for plotting the dataset(s) (default = 'brgkywrgbkywrgbkywrgbkyw')
%                       alternatively, colors can be specified as nx3 matrix of rgb values
%   cfg.directionality = '', 'inflow' or 'outflow' specifies for
%                       connectivity measures whether the inflow into a
%                       node, or the outflow from a node is plotted. The
%                       (default) behavior of this option depends on the dimor
%                       of the input data (see below).
%
% For the plotting of directional connectivity data the cfg.directionality
% option determines what is plotted. The default value and the supported
% functionality depend on the dimord of the input data. If the input data
% is of dimord 'chan_chan_XXX', the value of directionality determines
% whether, given the reference channel(s), the columns (inflow), or rows
% (outflow) are selected for plotting. In this situation the default is
% 'inflow'. Note that for undirected measures, inflow and outflow should
% give the same output. If the input data is of dimord 'chancmb_XXX', the
% value of directionality determines whether the rows in data.labelcmb are
% selected. With 'inflow' the rows are selected if the refchannel(s) occur in
% the right column, with 'outflow' the rows are selected if the
% refchannel(s) occur in the left column of the labelcmb-field. Default in
% this case is '', which means that all rows are selected in which the
% refchannel(s) occur. This is to robustly support linearly indexed
% undirected connectivity metrics. In the situation where undirected
% connectivity measures are linearly indexed, specifying 'inflow' or
% 'outflow' can result in unexpected behavior.
%
% to facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% if you specify this option the input data will be read from a *.mat
% file on disk. this mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SINGLEPLOTTFR, FT_MULTIPLOTER, FT_MULTIPLOTTFR, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% Undocumented local options:
% cfg.zlim/xparam (set to a specific frequency range or time range [zmax zmin] for an average
% over the frequency/time bins for TFR data.  Use in conjunction with e.g. xparam = 'time', and cfg.parameter = 'powspctrm').
% cfg.preproc

% Copyright (C) 2003-2006, Ole Jensen
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
cfg = ft_checkconfig(cfg, 'unused',     {'cohtargetchannel'});
cfg = ft_checkconfig(cfg, 'renamedval', {'zlim', 'absmax', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedforward', 'outflow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'directionality', 'feedback',    'inflow'});
cfg = ft_checkconfig(cfg, 'renamed',    {'matrixside',     'directionality'});
cfg = ft_checkconfig(cfg, 'renamed',    {'channelindex',   'channel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'channelname',    'channel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'cohrefchannel',  'refchannel'});
cfg = ft_checkconfig(cfg, 'renamed',	  {'zparam',         'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated', {'xparam'});

% set the defaults
cfg.baseline        = ft_getopt(cfg, 'baseline',    'no');
cfg.trials          = ft_getopt(cfg, 'trials',      'all', 1);
cfg.xlim            = ft_getopt(cfg, 'xlim',        'maxmin');
cfg.ylim            = ft_getopt(cfg, 'ylim',        'maxmin');
cfg.zlim            = ft_getopt(cfg, 'zlim',        'maxmin');
cfg.comment         = ft_getopt(cfg, 'comment',     strcat([date '\n']));
cfg.axes            = ft_getopt(cfg,' axes',        'yes');
cfg.fontsize        = ft_getopt(cfg, 'fontsize',    8);
cfg.graphcolor      = ft_getopt(cfg, 'graphcolor',  'brgkywrgbkywrgbkywrgbkyw');
cfg.hotkeys         = ft_getopt(cfg, 'hotkeys', 'no');
cfg.interactive     = ft_getopt(cfg, 'interactive',  'yes');
cfg.renderer        = ft_getopt(cfg, 'renderer',     []);
cfg.maskparameter   = ft_getopt(cfg, 'maskparameter',[]);
cfg.linestyle       = ft_getopt(cfg, 'linestyle',    '-');
cfg.linewidth       = ft_getopt(cfg, 'linewidth',    0.5);
cfg.maskstyle       = ft_getopt(cfg, 'maskstyle',    'box');
cfg.channel         = ft_getopt(cfg, 'channel',      'all');
cfg.directionality  = ft_getopt(cfg, 'directionality',   []);
cfg.figurename      = ft_getopt(cfg, 'figurename',       []);
cfg.preproc         = ft_getopt(cfg, 'preproc', []);
cfg.frequency       = ft_getopt(cfg, 'frequency', 'all'); % needed for frequency selection with TFR data
cfg.latency         = ft_getopt(cfg, 'latency', 'all'); % needed for latency selection with TFR data, FIXME, probably not used


Ndata = numel(varargin);

% interactive plotting is not allowed with more than 1 input
% if Ndata >1 && strcmp(cfg.interactive, 'yes')
%   error('interactive plotting is not supported with more than 1 input data set');
% end

% FIXME rename directionality and cohrefchannel in more meaningful options
if ischar(cfg.graphcolor)
  graphcolor = ['k' cfg.graphcolor];
elseif isnumeric(cfg.graphcolor)
  graphcolor = [0 0 0; cfg.graphcolor];
end

% check for linestyle being a cell-array, check it's length, and lengthen it if does not have enough styles in it
if ischar(cfg.linestyle)
  cfg.linestyle = {cfg.linestyle};
end

if Ndata  > 1
  if (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) > 1)
    error('either specify cfg.linestyle as a cell-array with one cell for each dataset, or only specify one linestyle')
  elseif (length(cfg.linestyle) < Ndata ) && (length(cfg.linestyle) == 1)
    tmpstyle = cfg.linestyle{1};
    cfg.linestyle = cell(Ndata , 1);
    for idataset = 1:Ndata
      cfg.linestyle{idataset} = tmpstyle;
    end
  end
end

% ensure that the input is correct, also backward compatibility with old data structures:
dtype = cell(Ndata, 1);
for i=1:Ndata
  % check if the input data is valid for this function
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'timelock', 'freq'});
  dtype{i}    = ft_datatype(varargin{i});

  % this is needed for correct treatment of graphcolor later on
  if nargin>1,
    if ~isempty(inputname(i+1))
      iname{i+1} = inputname(i+1);
    else
      iname{i+1} = ['input',num2str(i,'%02d')];
    end
  else
    iname{i+1} = cfg.inputfile{i};
  end
end

if Ndata >1,
  if ~all(strcmp(dtype{1}, dtype))
    error('input data are of different type; this is not supported');
  end
end
dtype  = dtype{1};
dimord = varargin{1}.dimord;
dimtok = tokenize(dimord, '_');

% ensure that the preproc specific options are located in the cfg.preproc
% substructure, but also ensure that the field 'refchannel' is present at the
% highest level in the structure. This is a little hack by JM because the field
% refchannel can also refer to the plotting of a connectivity metric. Also,
% the freq2raw conversion does not work at all in the call to ft_preprocessing.
% Therefore, for now, the preprocessing will not be done when there is freq
% data in the input. A more generic solution should be considered.

if isfield(cfg, 'refchannel'), refchannelincfg = cfg.refchannel; end
if ~any(strcmp({'freq','freqmvar'},dtype)),
  cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
end
if exist('refchannelincfg', 'var'), cfg.refchannel  = refchannelincfg; end

if ~isempty(cfg.preproc)
  % preprocess the data, i.e. apply filtering, baselinecorrection, etc.
  fprintf('applying preprocessing options\n');
  if ~isfield(cfg.preproc, 'feedback')
    cfg.preproc.feedback = cfg.interactive;
  end
  for i=1:Ndata
    varargin{i} = ft_preprocessing(cfg.preproc, varargin{i});
  end
end

% set x/y/parameter defaults according to datatype and dimord
switch dtype
  case 'timelock'
    xparam = 'time';
    yparam = '';
    cfg.parameter = ft_getopt(cfg,  'parameter', 'avg');
  case 'freq'
    if sum(ismember(dimtok, 'time'))
      xparam = 'time';
      yparam = 'freq';
      cfg.parameter = ft_getopt(cfg,  'parameter', 'powspctrm');
    elseif sum(ismember(dimtok, 'time'))
      xparam = 'freq';
      yparam = 'time';
      cfg.parameter = ft_getopt(cfg,  'parameter', 'powspctrm');
    else
      xparam = 'freq';
      yparam = '';
      cfg.parameter = ft_getopt(cfg,  'parameter', 'powspctrm');
    end
  case 'comp'
    % not supported
  otherwise
    % not supported
end

% user specified own fields, but no yparam (which is not asked in help)
if exist('xparam', 'var') && isfield(cfg, 'parameter') && ~exist('yparam', 'var')
  yparam = '';
end

if isfield(cfg, 'channel') && isfield(varargin{1}, 'label')
  cfg.channel = ft_channelselection(cfg.channel, varargin{1}.label);
elseif isfield(cfg, 'channel') && isfield(varargin{1}, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(varargin{1}.labelcmb(:)));
end

% check whether rpt/subj is present and remove if necessary and whether
hasrpt = sum(ismember(dimtok, {'rpt' 'subj'}));
if strcmp(dtype, 'timelock') && hasrpt,
  tmpcfg        = [];
  tmpcfg.trials = cfg.trials;
  for i=1:Ndata
    varargin{i} = ft_timelockanalysis(tmpcfg, varargin{i});
  end
  if ~strcmp(cfg.parameter, 'avg')
    % rename avg back into the parameter
    varargin{i}.(cfg.parameter) = varargin{i}.avg;
    varargin{i}                 = rmfield(varargin{i}, 'avg');
  end
  dimord        = varargin{1}.dimord;
  dimtok        = tokenize(dimord, '_');
elseif strcmp(dtype, 'freq') && hasrpt,
  % this also deals with fourier-spectra in the input
  % or with multiple subjects in a frequency domain stat-structure
  % on the fly computation of coherence spectrum is not supported
  for i=1:Ndata
    if isfield(varargin{i}, 'crsspctrm'),
      varargin{i} = rmfield(varargin{i}, 'crsspctrm');
    end
  end

  tmpcfg           = [];
  tmpcfg.trials    = cfg.trials;
  tmpcfg.jackknife = 'no';
  for i=1:Ndata
    if isfield(cfg, 'parameter') && ~strcmp(cfg.parameter,'powspctrm')
      % freqdesctiptives will only work on the powspctrm field
      % hence a temporary copy of the data is needed
      tempdata.dimord    = varargin{i}.dimord;
      tempdata.freq      = varargin{i}.freq;
      tempdata.label     = varargin{i}.label;
      tempdata.powspctrm = varargin{i}.(cfg.parameter);
      if isfield(varargin{i}, 'cfg') tempdata.cfg = varargin{i}.cfg; end
      tempdata           = ft_freqdescriptives(tmpcfg, tempdata);
      varargin{i}.(cfg.parameter)  = tempdata.powspctrm;
      clear tempdata
    else
      varargin{i} = ft_freqdescriptives(tmpcfg, varargin{i});
    end
  end
  dimord = varargin{1}.dimord;
  dimtok = tokenize(dimord, '_');
end

% apply baseline correction
if ~strcmp(cfg.baseline, 'no')
  for i=1:Ndata
    if strcmp(dtype, 'timelock') && strcmp(xparam, 'time')
      varargin{i} = ft_timelockbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'time')
      varargin{i} = ft_freqbaseline(cfg, varargin{i});
    elseif strcmp(dtype, 'freq') && strcmp(xparam, 'freq')
      error('baseline correction is not supported for spectra without a time dimension');
    else
      warning('baseline correction not applied, please set xparam');
    end
  end
end

% handle the bivariate case

% check for bivariate metric with 'chan_chan' in the dimord
selchan = strmatch('chan', dimtok);
isfull  = length(selchan)>1;

% check for bivariate metric with a labelcmb
haslabelcmb = isfield(varargin{1}, 'labelcmb');

if (isfull || haslabelcmb) && (isfield(varargin{1}, cfg.parameter) && ~strcmp(cfg.parameter, 'powspctrm'))
  % a reference channel is required:
  if ~isfield(cfg, 'refchannel')
    error('no reference channel is specified');
  end

  % check for refchannel being part of selection
  if ~strcmp(cfg.refchannel,'gui')
    if haslabelcmb
      cfg.refchannel = ft_channelselection(cfg.refchannel, unique(varargin{1}.labelcmb(:)));
    else
      cfg.refchannel = ft_channelselection(cfg.refchannel, varargin{1}.label);
    end
    if (isfull      && ~any(ismember(varargin{1}.label, cfg.refchannel))) || ...
        (haslabelcmb && ~any(ismember(varargin{1}.labelcmb(:), cfg.refchannel)))
      error('cfg.refchannel is a not present in the (selected) channels)')
    end
  end

  % interactively select the reference channel
  if strcmp(cfg.refchannel, 'gui')
    error('cfg.refchannel = ''gui'' is not supported in ft_singleplotER');
  end

  for i=1:Ndata
    if ~isfull,
      % convert 2-dimensional channel matrix to a single dimension:
      if isempty(cfg.directionality)
        sel1 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:,2)));
        sel2 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:,1)));
      elseif strcmp(cfg.directionality, 'outflow')
        sel1 = [];
        sel2 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:,1)));
      elseif strcmp(cfg.directionality, 'inflow')
        sel1 = find(strcmp(cfg.refchannel, varargin{i}.labelcmb(:,2)));
        sel2 = [];
      end
      fprintf('selected %d channels for %s\n', length(sel1)+length(sel2), cfg.parameter);
      if length(sel1)+length(sel2)==0
        error('there are no channels selected for plotting: you may need to look at the specification of cfg.directionality');
      end
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)([sel1;sel2],:,:);
      varargin{i}.label     = [varargin{i}.labelcmb(sel1,1);varargin{i}.labelcmb(sel2,2)];
      varargin{i}.labelcmb  = varargin{i}.labelcmb([sel1;sel2],:);
      varargin{i}           = rmfield(varargin{i}, 'labelcmb');
    else
      % general case
      sel               = match_str(varargin{i}.label, cfg.refchannel);
      siz               = [size(varargin{i}.(cfg.parameter)) 1];
      if strcmp(cfg.directionality, 'inflow') || isempty(cfg.directionality)
        %the interpretation of 'inflow' and 'outflow' depend on
        %the definition in the bivariate representation of the data
        %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(:,sel,:),2),[siz(1) 1 siz(3:end)]);
        sel1 = 1:siz(1);
        sel2 = sel;
        meandir = 2;
      elseif strcmp(cfg.directionality, 'outflow')
        %data.(cfg.parameter) = reshape(mean(data.(cfg.parameter)(sel,:,:),1),[siz(1) 1 siz(3:end)]);
        sel1 = sel;
        sel2 = 1:siz(1);
        meandir = 1;

      elseif strcmp(cfg.directionality, 'ff-fd')
        error('cfg.directionality = ''ff-fd'' is not supported anymore, you have to manually subtract the two before the call to ft_singleplotER');
      elseif strcmp(cfg.directionality, 'fd-ff')
        error('cfg.directionality = ''fd-ff'' is not supported anymore, you have to manually subtract the two before the call to ft_singleplotER');
      end %if directionality
    end %if ~isfull
  end %for i
end %handle the bivariate data

% get physical min/max range of x
if strcmp(cfg.xlim,'maxmin')
  % find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:Ndata
    xmin = min([xmin varargin{i}.(xparam)]);
    xmax = max([xmax varargin{i}.(xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% get the index of the nearest bin
for i=1:Ndata
  xidmin(i,1) = nearest(varargin{i}.(xparam), xmin);
  xidmax(i,1) = nearest(varargin{i}.(xparam), xmax);
end

if strcmp('freq', yparam) && strcmp('freq', dtype)
  tmpcfg = keepfields(cfg, {'parameter'});
  tmpcfg.avgoverfreq = 'yes';
  tmpcfg.frequency   = cfg.frequency;%cfg.zlim;
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
elseif strcmp('time', yparam) && strcmp('freq', dtype)
  tmpcfg = keepfields(cfg, {'parameter'});
  tmpcfg.avgovertime = 'yes';
  tmpcfg.latency     = cf.latency;%cfg.zlim;
  [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
  % restore the provenance information
  [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
end

cla
hold on;
colorlabels = [];

% plot each data set:
for i=1:Ndata
  if isfield(varargin{1}, 'label')
    selchannel = ft_channelselection(cfg.channel, varargin{i}.label);
  elseif isfield(varargin{1}, 'labelcmb')
    selchannel = ft_channelselection(cfg.channel, unique(varargin{i}.labelcmb(:)));
  else
    error('the input data does not contain a label or labelcmb-field');
  end

  % make vector dat with one value for each channel
  dat  = varargin{i}.(cfg.parameter);
  % get dimord dimensions
  dims = textscan(varargin{i}.dimord,'%s', 'Delimiter', '_');
  dims = dims{1};
  ydim = find(strcmp(yparam, dims));
  xdim = find(strcmp(xparam, dims));
  zdim = setdiff(1:ndims(dat), [ydim xdim]);
  % and permute to make sure that dimensions are in the correct order
  dat = permute(dat, [zdim(:)' ydim xdim]);


  xval = varargin{i}.(xparam);

  % take subselection of channels
  % this works for bivariate data with labelcmb because at this point the
  % data has a label-field
  sellab = match_str(varargin{i}.label, selchannel);

  %     if ~isempty(yparam)
  %         if isfull
  %             dat = dat(sel1, sel2, ymin:ymax, xidmin(i):xidmax(i));
  %             dat = nanmean(nanmean(dat, meandir), 3);
  %             siz = size(dat);
  %             %fixmedat = reshape(dat, [siz(1:2) siz(4)]);
  %             dat = reshape(dat, [siz(1) siz(3)]);
  %             dat = dat(sellab, :);
  %         elseif haslabelcmb
  %             dat = dat(sellab, ymin:ymax, xidmin(i):xidmax(i));
  %             dat = nanmean(dat, 2);
  %             siz = size(dat);
  %             dat = reshape(dat, [siz(1) siz(3)]);
  %         else
  %             dat = dat(sellab, ymin:ymax, xidmin(i):xidmax(i));
  %             dat = nanmean(nanmean(dat, 3), 2);
  %             siz = size(dat);
  %             dat = reshape(dat, [siz(1) siz(3)]);
  %         end
  %     else
  if isfull
    dat = dat(sel1, sel2, xidmin(i):xidmax(i));
    dat = nanmean(dat, meandir);
    siz = size(dat);
    siz(find(siz(1:2)==1)) = [];
    dat = reshape(dat, siz);
    dat = dat(sellab, :);
  elseif haslabelcmb
    dat = dat(sellab, xidmin(i):xidmax(i));
  else
    dat = dat(sellab, xidmin(i):xidmax(i));
  end
  %     end
  xval       = xval(xidmin(i):xidmax(i));
  datavector = reshape(mean(dat, 1), [1 numel(xval)]); % average over channels

  % make mask
  if ~isempty(cfg.maskparameter)
    datmask = varargin{i}.(cfg.maskparameter)(sellab,:);
    if size(datmask,2)>1
      datmask = datmask(:,xidmin(i):xidmax(i));
    else
      datmask = datmask(xidmin(i):xidmax(i));
    end
    maskdatavector = reshape(mean(datmask,1), [1 numel(xval)]);
  else
    maskdatavector = [];
  end

  if Ndata  > 1
    if ischar(graphcolor);        colorlabels = [colorlabels iname{i+1} '=' graphcolor(i+1) '\n'];
    elseif isnumeric(graphcolor); colorlabels = [colorlabels iname{i+1} '=' num2str(graphcolor(i+1,:)) '\n'];
    end
  end

  if ischar(graphcolor);        color = graphcolor(i+1);
  elseif isnumeric(graphcolor); color = graphcolor(i+1,:);
  end

  % update ymin and ymax for the current data set:
  if ischar(cfg.ylim)
    if i==1
      ymin = [];
      ymax = [];
    end
    if strcmp(cfg.ylim,'maxmin')
      % select the channels in the data that match with the layout:
      ymin = min([ymin min(datavector)]);
      ymax = max([ymax max(datavector)]);
    elseif strcmp(cfg.ylim,'maxabs')
      ymax = max([ymax max(abs(datavector))]);
      ymin = -ymax;
    elseif strcmp(cfg.ylim,'zeromax')
      ymin = 0;
      ymax = max([ymax max(datavector)]);
    elseif strcmp(cfg.ylim,'minzero')
      ymin = min([ymin min(datavector)]);
      ymax = 0;
    end;
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end


  % only plot the mask once, for the first line (it's the same anyway for
  % all lines, and if plotted multiple times, it will overlay the others
  if i>1 && strcmp(cfg.maskstyle, 'box')
    ft_plot_vector(xval, datavector, 'style', cfg.linestyle{i}, 'color', color, ...
      'linewidth', cfg.linewidth, 'hlim', cfg.xlim, 'vlim', cfg.ylim);
  else
    ft_plot_vector(xval, datavector, 'style', cfg.linestyle{i}, 'color', color, ...
      'highlight', maskdatavector, 'highlightstyle', cfg.maskstyle, 'linewidth', cfg.linewidth, ...
      'hlim', cfg.xlim, 'vlim', cfg.ylim);
  end
end

% set xlim and ylim:
xlim([xmin xmax]);
ylim([ymin ymax]);

% adjust mask box extents to ymin/ymax
if ~isempty(cfg.maskparameter)
  ptchs = findobj(gcf,'type','patch');
  for i = 1:length(ptchs)
    YData = get(ptchs(i),'YData');
    YData(YData == min(YData)) = ymin;
    YData(YData == max(YData)) = ymax;
    set(ptchs(i),'YData',YData);
  end
end

if strcmp('yes',cfg.hotkeys)
  %  attach data and cfg to figure and attach a key listener to the figure
  set(gcf, 'keypressfcn', {@key_sub, xmin, xmax, ymin, ymax})
end

if isfield(cfg, 'dataname')
  dataname = cfg.dataname;
elseif nargin > 1
  dataname = inputname(2);
  cfg.dataname = {inputname(2)};
  for k = 2:Ndata
    dataname = [dataname ', ' inputname(k+1)];
    cfg.dataname{end+1} = inputname(k+1);
  end
else
  dataname = cfg.inputfile;
end

% set the figure window title, add the channel labels if number is small
if isempty(get(gcf,'Name'))
  if length(sellab) < 5
    chans = join_str(',', cfg.channel);
  else
    chans = '<multiple channels>';
  end
  if isempty(cfg.figurename)
    set(gcf, 'Name', sprintf('%d: %s: %s (%s)', double(gcf), mfilename, join_str(', ',dataname), chans));
    set(gcf, 'NumberTitle', 'off');
  else
    set(gcf, 'name', cfg.figurename);
    set(gcf, 'NumberTitle', 'off');
  end
end

% make the figure interactive
if strcmp(cfg.interactive, 'yes')
  % add the dataname to the figure
  % this is used in the callbacks
  info          = guidata(gcf);
  info.dataname = dataname;
  guidata(gcf, info);
  % attach data to the figure with the current axis handle as a name
  dataname = fixname(num2str(double(gca)));
  setappdata(gcf,dataname,varargin);
  set(gcf, 'windowbuttonupfcn',     {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg}, 'event', 'windowbuttonupfcn'});
  set(gcf, 'windowbuttondownfcn',   {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg}, 'event', 'windowbuttondownfcn'});
  set(gcf, 'windowbuttonmotionfcn', {@ft_select_range, 'multiple', false, 'yrange', false, 'callback', {@select_topoplotER, cfg}, 'event', 'windowbuttonmotionfcn'});
end

% create title text containing channel name(s) and channel number(s):
if length(sellab) == 1
  t = [char(cfg.channel) ' / ' num2str(sellab) ];
else
  t = sprintf('mean(%0s)', join_str(',', cfg.channel));
end
h = title(t,'fontsize', cfg.fontsize);

% set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

if false
  % FIXME this is for testing purposes
  % Define a context menu; it is not attached to anything
  cmlines = uicontextmenu;
  % Define the context menu items and install their callbacks
  uimenu(cmlines, 'Label', 'dashed', 'Callback', 'set(gco, ''LineStyle'', ''--'')');
  uimenu(cmlines, 'Label', 'dotted', 'Callback', 'set(gco, ''LineStyle'', '':'')');
  uimenu(cmlines, 'Label', 'solid',  'Callback', 'set(gco, ''LineStyle'', ''-'')');
  % Locate line objects
  hlines = findall(gca, 'Type', 'line');
  % Attach the context menu to each line
  for line = 1:length(hlines)
    set(hlines(line), 'uicontextmenu', cmlines)
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance

% add a menu to the figure, but only if the current figure does not have subplots
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
delete(findobj(gcf, 'type', 'uimenu', 'label', 'FieldTrip'));
if numel(findobj(gcf, 'type', 'axes', '-not', 'tag', 'ft-colorbar')) <= 1
  ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
  uimenu(ftmenu, 'Label', 'Show pipeline',  'Callback', {@menu_pipeline, cfg});
  uimenu(ftmenu, 'Label', 'About',  'Callback', @menu_about);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting a time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_topoplotER(cfg, varargin)
% first to last callback-input of ft_select_range is range
% last callback-input of ft_select_range is contextmenu label, if used
range = varargin{end-1};
varargin = varargin(1:end-2); % remove range and last

% get appdata belonging to current axis
dataname = fixname(num2str(double(gca)));
data = getappdata(gcf, dataname);

if isfield(cfg, 'inputfile')
  % the reading has already been done and varargin contains the data
  cfg = rmfield(cfg, 'inputfile');
end
if isfield(cfg, 'showlabels')
  % this is not allowed in topoplotER
  cfg = rmfield(cfg, 'showlabels');
end
% make sure the topo displays all channels, not just the ones in this singleplot
cfg.channel = 'all';
cfg.comment = 'auto';
cfg.xlim = range(1:2);
% put data name in here, this cannot be resolved by other means
info = guidata(gcf);
cfg.dataname = info.dataname;
% if user specified a ylim, copy it over to the zlim of topoplot
if isfield(cfg, 'ylim')
  cfg.zlim = cfg.ylim;
  cfg = rmfield(cfg, 'ylim');
end
fprintf('selected cfg.xlim = [%f %f]\n', cfg.xlim(1), cfg.xlim(2));
p = get(gcf, 'position');
f = figure;
set(f, 'position', p);
ft_topoplotER(cfg, data{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata, varargin)
xlimits = xlim;
ylimits = ylim;
incr_x = abs(xlimits(2) - xlimits(1)) /10;
incr_y = abs(ylimits(2) - ylimits(1)) /10;

% TRANSLATE by 10%
if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'leftarrow')
  xlim([xlimits(1)+incr_x xlimits(2)+incr_x])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'rightarrow')
  xlim([xlimits(1)-incr_x xlimits(2)-incr_x])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'uparrow')
  ylim([ylimits(1)-incr_y ylimits(2)-incr_y])
elseif length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{:},'control') && strcmp(eventdata.Key,'downarrow')
  ylim([ylimits(1)+incr_y ylimits(2)+incr_y])
  % ZOOM by 10%
elseif strcmp(eventdata.Key,'leftarrow')
  xlim([xlimits(1)-incr_x xlimits(2)+incr_x])
elseif strcmp(eventdata.Key,'rightarrow')
  xlim([xlimits(1)+incr_x xlimits(2)-incr_x])
elseif strcmp(eventdata.Key,'uparrow')
  ylim([ylimits(1)-incr_y ylimits(2)+incr_y])
elseif strcmp(eventdata.Key,'downarrow')
  ylim([ylimits(1)+incr_y ylimits(2)-incr_y])
  % resort to minmax of data for x-axis and y-axis
elseif strcmp(eventdata.Key,'m')
  xlim([varargin{1} varargin{2}])
  ylim([varargin{3} varargin{4}])
end
