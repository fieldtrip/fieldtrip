function [cfg] = ft_clusterplot(cfg, stat)

% FT_CLUSTERPLOT plots a series of topographies with highlighted clusters.
%
% Use as
%   ft_clusterplot(cfg, stat)
% where the input data is obtained from FT_TIMELOCKSTATISTICS or FT_FREQSTATISTICS.
%
% The configuration options can be
%   cfg.alpha                     = number, highest cluster p-value to be plotted max 0.3 (default = 0.05)
%   cfg.highlightseries           = 1x5 cell-array, highlight option series  with 'on', 'labels' or 'numbers' (default {'on', 'on', 'on', 'on', 'on'} for p < [0.01 0.05 0.1 0.2 0.3]
%   cfg.highlightsymbolseries     = 1x5 vector, highlight marker symbol series (default ['*', 'x', '+', 'o', '.'] for p < [0.01 0.05 0.1 0.2 0.3]
%   cfg.highlightsizeseries       = 1x5 vector, highlight marker size series   (default [6 6 6 6 6] for p < [0.01 0.05 0.1 0.2 0.3])
%   cfg.highlightcolorpos         = color of highlight marker for positive clusters (default = [0 0 0])
%   cfg.highlightcolorneg         = color of highlight marker for negative clusters (default = [0 0 0])
%   cfg.subplotsize               = layout of subplots ([h w], default [3 5])
%   cfg.saveaspng                 = string, filename of the output figures (default = 'no')
%   cfg.visible                   = string, 'on' or 'off' whether figure will be visible (default = 'on')
%
% You can also specify all cfg options that apply to FT_TOPOPLOTER or FT_TOPOPLOTTFR,
% except for cfg.xlim, any of the highlight options, cfg.comment and cfg.commentpos.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_TOPOPLOTTFR, FT_TOPOPLOTER, FT_MOVIEPLOTTFR, FT_MOVIEPLOTER

% Copyright (C) 2007, F.C. Donders Centre, Ingrid Nieuwenhuis
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
ft_preamble loadvar stat
ft_preamble provenance stat
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

ws = ft_warning('off', 'FieldTrip:getdimord:warning_dimord_could_not_be_determined');

% check if the input data is valid for this function
stat = ft_checkdata(stat, 'datatype', {'timelock', 'freq'}, 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'hlmarkerseries',       'highlightsymbolseries'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlmarkersizeseries',   'highlightsizeseries'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlcolorpos',           'highlightcolorpos'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hlcolorneg',           'highlightcolorneg'});
cfg = ft_checkconfig(cfg, 'renamed',     {'zparam', 'parameter'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'hllinewidthseries'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'xparam', 'yparam'});

% added several forbidden options
cfg = ft_checkconfig(cfg, 'forbidden',  {'highlight', ...
  'highlightchannel', ...
  'highlightsymbol', ...
  'highlightcolor', ...
  'highlightsize', ...
  'highlightfontsize', ...
  'xlim', ...
  'comment', ...
  'commentpos'});

% set the defaults
cfg.highlightseries         = ft_getopt(cfg, 'highlightseries',         {'on', 'on', 'on', 'on', 'on'});
cfg.highlightsymbolseries   = ft_getopt(cfg, 'highlightsymbolseries',   ['*', 'x', '+', 'o', '.']);
cfg.highlightsizeseries     = ft_getopt(cfg, 'highlightsizeseries',     [6 6 6 6 6]);
cfg.hllinewidthseries       = ft_getopt(cfg, 'hllinewidthseries',       [1 1 1 1 1]);
cfg.highlightfontsizeseries = ft_getopt(cfg, 'highlightfontsizeseries', [8 8 8 8 8]);
cfg.highlightcolorpos       = ft_getopt(cfg, 'highlightcolorpos',       [0 0 0]);
cfg.highlightcolorneg       = ft_getopt(cfg, 'highlightcolorneg',       [0 0 0]);
cfg.marker                  = ft_getopt(cfg, 'marker',                  'off');
cfg.alpha                   = ft_getopt(cfg, 'alpha',                   0.05);
cfg.parameter               = ft_getopt(cfg, 'parameter',               'stat');
cfg.saveaspng               = ft_getopt(cfg, 'saveaspng',               'no');
cfg.subplotsize             = ft_getopt(cfg, 'subplotsize',             [3 5]);
cfg.feedback                = ft_getopt(cfg, 'feedback',                'text');
cfg.visible                 = ft_getopt(cfg, 'visible',                 'on');

% error if cfg.highlightseries is not a cell, for possible confusion with cfg-options
if ~iscell(cfg.highlightseries)
  ft_error('cfg.highlightseries should be a cell-array of strings')
end

% get the options that are specific for topoplotting
cfgtopo = keepfields(cfg, {'parameter', 'marker', 'markersymbol', 'markercolor', 'markersize', 'markerfontsize', 'style', 'gridscale', 'interplimits', 'interpolation', 'contournum', 'colorbar', 'shading', 'zlim'});
% prepare the layout, this only has to be done once
cfgtopo.layout = ft_prepare_layout(cfg, stat);
cfgtopo.showcallinfo = 'no';
cfgtopo.feedback = 'no';

% handle with the data, it should be 1D or 2D
dimord = getdimord(stat, cfg.parameter);
dimtok = tokenize(dimord, '_');
dimsiz = getdimsiz(stat, cfg.parameter);
dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions

switch dimord
  case 'chan'
    is2D = false;
    
  case 'chan_time'
    is2D = true;
    
  case 'chan_freq'
    is2D = true;
    
  case 'chan_freq_time'
    % no more than two dimensions are supported, we can ignore singleton dimensions
    is2D = true;
    if dimsiz(2)==1
      stat = rmfield(stat, 'freq');
      stat.dimord = 'chan_time';
      % remove the singleton dimension in the middle
      stat.(cfg.parameter) = reshape(stat.(cfg.parameter),dimsiz([1 3]));
      if isfield(stat, 'posclusterslabelmat')
        stat.posclusterslabelmat = reshape(stat.posclusterslabelmat, dimsiz([1 3]));
      end
      if isfield(stat, 'negclusterslabelmat')
        stat.negclusterslabelmat = reshape(stat.negclusterslabelmat, dimsiz([1 3]));
      end
    elseif dimsiz(3)==1
      stat = rmfield(stat, 'time');
      stat.dimord = 'chan_freq';
      % no need to remove the singleton dimension at the end
    else
      ft_error('this only works if either frequency or time is a singleton dimension');
    end
    
  otherwise
    ft_error('unsupported dimord %s', dimord);
end % switch dimord

% these are not valid any more
clear dimord dimsiz

% this determines the labels in the figure
hastime = isfield(stat, 'time');
hasfreq = isfield(stat, 'freq');

% use the vector time, even though the 2nd dimension might be freq
if hastime
  time = stat.time;
elseif hasfreq
  time = stat.freq;
end

if issubfield(stat, 'cfg.correcttail') && ((strcmp(stat.cfg.correcttail, 'alpha') || strcmp(stat.cfg.correcttail, 'prob')) && (stat.cfg.tail == 0));
  if ~(cfg.alpha >= stat.cfg.alpha)
    ft_warning(['the pvalue you plot: cfg.alpha = ' num2str(cfg.alpha) ' is higher than the correcttail option you tested: stat.cfg.alpha = ' num2str(stat.cfg.alpha)]);
  end
end

% find significant clusters
sigpos = [];
signeg = [];
haspos = isfield(stat, 'posclusters');
hasneg = isfield(stat, 'negclusters');

if haspos == 0 && hasneg == 0
  fprintf('%s\n', 'no significant clusters in data; nothing to plot')
else
  if haspos
    for iPos = 1:length(stat.posclusters)
      sigpos(iPos) = stat.posclusters(iPos).prob < cfg.alpha;
    end
  end
  if hasneg
    for iNeg = 1:length(stat.negclusters)
      signeg(iNeg) = stat.negclusters(iNeg).prob < cfg.alpha;
    end
  end
  sigpos  = find(sigpos == 1);
  signeg  = find(signeg == 1);
  Nsigpos = length(sigpos);
  Nsigneg = length(signeg);
  Nsigall = Nsigpos + Nsigneg;
  
  if Nsigall == 0
    ft_error('no clusters present with a p-value lower than the specified alpha, nothing to plot')
  end
  
  % make clusterslabel matrix per significant cluster
  if haspos
    posCLM = stat.posclusterslabelmat;
    sigposCLM = zeros(size(posCLM));
    probpos = [];
    for iPos = 1:length(sigpos)
      sigposCLM(:,:,iPos) = (posCLM == sigpos(iPos));
      probpos(iPos) = stat.posclusters(iPos).prob;
      hlsignpos(iPos) = prob2hlsign(probpos(iPos), cfg.highlightsymbolseries);
    end
  else
    posCLM = [];
    sigposCLM = [];
    probpos = [];
  end
  
  if hasneg
    negCLM = stat.negclusterslabelmat;
    signegCLM = zeros(size(negCLM));
    probneg = [];
    for iNeg = 1:length(signeg)
      signegCLM(:,:,iNeg) = (negCLM == signeg(iNeg));
      probneg(iNeg) = stat.negclusters(iNeg).prob;
      hlsignneg(iNeg) = prob2hlsign(probneg(iNeg), cfg.highlightsymbolseries);
    end
  else % no negative clusters
    negCLM = [];
    signegCLM = [];
    probneg = [];
  end
  
  fprintf('There are %d clusters smaller than alpha (%g)\n', Nsigall, cfg.alpha);
  
  if is2D
    % define time or freq window per cluster
    for iPos = 1:length(sigpos)
      possum_perclus = sum(sigposCLM(:,:,iPos),1); %sum over chans for each time- or freq-point
      ind_min = find(possum_perclus~=0, 1 );
      ind_max = find(possum_perclus~=0, 1, 'last' );
      time_perclus = [time(ind_min) time(ind_max)];
      if hastime
        fprintf('%s%s%s%s%s%s%s%s%s%s%s\n', 'Positive cluster: ',num2str(sigpos(iPos)), ', pvalue: ',num2str(probpos(iPos)), ' (',hlsignpos(iPos), ')', ', t = ',num2str(time_perclus(1)), ' to ',num2str(time_perclus(2)))
      elseif hasfreq
        fprintf('%s%s%s%s%s%s%s%s%s%s%s\n', 'Positive cluster: ',num2str(sigpos(iPos)), ', pvalue: ',num2str(probpos(iPos)), ' (',hlsignpos(iPos), ')', ', f = ',num2str(time_perclus(1)), ' to ',num2str(time_perclus(2)))
      end
    end
    for iNeg = 1:length(signeg)
      negsum_perclus = sum(signegCLM(:,:,iNeg),1);
      ind_min = find(negsum_perclus~=0, 1 );
      ind_max = find(negsum_perclus~=0, 1, 'last' );
      time_perclus = [time(ind_min) time(ind_max)];
      if hastime
        time_perclus = [time(ind_min) time(ind_max)];
        fprintf('%s%s%s%s%s%s%s%s%s%s%s\n', 'Negative cluster: ',num2str(signeg(iNeg)), ', pvalue: ',num2str(probneg(iNeg)), ' (',hlsignneg(iNeg), ')', ', t = ',num2str(time_perclus(1)), ' to ',num2str(time_perclus(2)))
      elseif hasfreq
        fprintf('%s%s%s%s%s%s%s%s%s%s%s\n', 'Negative cluster: ',num2str(signeg(iNeg)), ', pvalue: ',num2str(probneg(iNeg)), ' (',hlsignneg(iNeg), ')', ', f = ',num2str(time_perclus(1)), ' to ',num2str(time_perclus(2)))
      end
    end
    
    % define time- or freq-window containing all significant clusters
    possum = sum(sigposCLM,3); %sum over Chans for timevector
    possum = sum(possum,1);
    negsum = sum(signegCLM,3);
    negsum = sum(negsum,1);
    
    if haspos && hasneg
      allsum = possum + negsum;
    elseif haspos
      allsum = possum;
    else
      allsum = negsum;
    end
    
    ind_timewin_min = find(allsum~=0, 1 );
    ind_timewin_max = find(allsum~=0, 1, 'last' );
    timewin = time(ind_timewin_min:ind_timewin_max);
    
  else
    for iPos = 1:length(sigpos)
      fprintf('%s%s%s%s%s%s%s\n', 'Positive cluster: ',num2str(sigpos(iPos)), ', pvalue: ',num2str(probpos(iPos)), ' (',hlsignpos(iPos), ')')
    end
    for iNeg = 1:length(signeg)
      fprintf('%s%s%s%s%s%s%s\n', 'Negative cluster: ',num2str(signeg(iNeg)), ', pvalue: ',num2str(probneg(iNeg)), ' (',hlsignneg(iNeg), ')')
    end
  end
  
  % setup highlight options for all clusters and make comment for 1D data
  compos = [];
  comneg = [];
  for iPos = 1:length(sigpos)
    if stat.posclusters(sigpos(iPos)).prob < 0.01
      cfgtopo.highlight{iPos}         = cfg.highlightseries{1};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(1);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(1);
      cfgtopo.highlightfontsize{iPos} = cfg.highlightfontsizeseries(1);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.05
      cfgtopo.highlight{iPos}         = cfg.highlightseries{2};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(2);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(2);
      cfgtopo.highlightfontsize{iPos} = cfg.highlightfontsizeseries(2);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.1
      cfgtopo.highlight{iPos}         = cfg.highlightseries{3};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(3);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(3);
      cfgtopo.highlightfontsize{iPos} = cfg.highlightfontsizeseries(3);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.2
      cfgtopo.highlight{iPos}         = cfg.highlightseries{4};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(4);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(4);
      cfgtopo.highlightfontsize{iPos} = cfg.highlightfontsizeseries(4);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.3
      cfgtopo.highlight{iPos}         = cfg.highlightseries{5};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(5);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(5);
      cfgtopo.highlightfontsize{iPos} = cfg.highlightfontsizeseries(5);
    end
    cfgtopo.highlightcolor{iPos}        = cfg.highlightcolorpos;
    compos = strcat(compos,cfgtopo.highlightsymbol{iPos}, 'p=',num2str(probpos(iPos)), ' '); % make comment, only used for 1D data
  end
  
  for iNeg = 1:length(signeg)
    if stat.negclusters(signeg(iNeg)).prob < 0.01
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{1};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(1);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(1);
      cfgtopo.highlightfontsize{length(sigpos)+iNeg} = cfg.highlightfontsizeseries(1);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.05
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{2};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(2);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(2);
      cfgtopo.highlightfontsize{length(sigpos)+iNeg} = cfg.highlightfontsizeseries(2);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.1
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{3};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(3);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(3);
      cfgtopo.highlightfontsize{length(sigpos)+iNeg} = cfg.highlightfontsizeseries(3);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.2
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{4};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(4);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(4);
      cfgtopo.highlightfontsize{length(sigpos)+iNeg} = cfg.highlightfontsizeseries(4);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.3
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{5};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(5);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(5);
      cfgtopo.highlightfontsize{length(sigpos)+iNeg} = cfg.highlightfontsizeseries(5);
    end
    cfgtopo.highlightcolor{length(sigpos)+iNeg}        = cfg.highlightcolorneg;
    comneg = strcat(comneg,cfgtopo.highlightsymbol{length(sigpos)+iNeg}, 'p=',num2str(probneg(iNeg)), ' '); % make comment, only used for 1D data
  end
  
  if is2D
    Npl = length(timewin);
  else
    Npl = 1;
  end
  
  numSubplots = prod(cfg.subplotsize);
  Nfig = ceil(Npl/numSubplots);
  
  % put channel indexes in list
  if is2D
    for iPl = 1:Npl
      for iPos = 1:length(sigpos)
        list{iPl}{iPos} = find(sigposCLM(:,ind_timewin_min+iPl-1,iPos) == 1);
      end
      for iNeg = 1:length(signeg)
        list{iPl}{length(sigpos)+iNeg} = find(signegCLM(:,ind_timewin_min+iPl-1,iNeg) == 1);
      end
    end
  else
    for iPl = 1:Npl
      for iPos = 1:length(sigpos)
        list{iPl}{iPos} = find(sigposCLM(:,iPos) == 1);
      end
      for iNeg = 1:length(signeg)
        list{iPl}{length(sigpos)+iNeg} = find(signegCLM(:,iNeg) == 1);
      end
    end
  end
  
  count = 0;
  ft_progress('init', cfg.feedback, 'making subplots...');
  ft_progress(count/Npl, 'making subplot %d from %d', count, Npl);
  
  % make plots
  for iPl = 1:Nfig
    f = figure('visible', cfg.visible);
    if is2D
      if iPl < Nfig
        for iT = 1:numSubplots
          PlN = (iPl-1)*numSubplots + iT; % plotnumber
          cfgtopo.xlim = [time(ind_timewin_min+PlN-1) time(ind_timewin_min+PlN-1)];
          cfgtopo.highlightchannel = list{PlN};
          if hastime
            cfgtopo.comment = strcat('time: ',num2str(time(ind_timewin_min+PlN-1)), ' s');
          elseif hasfreq
            cfgtopo.comment = strcat('freq: ',num2str(time(ind_timewin_min+PlN-1)), ' Hz');
          end
          cfgtopo.commentpos = 'title';
          subplot(cfg.subplotsize(1), cfg.subplotsize(2), iT);
          count = count+1;
          fprintf('making subplot %d from %d\n', count, Npl);
          ft_topoplotTFR(cfgtopo, stat);
        end
      elseif iPl == Nfig
        for iT = 1:Npl-(numSubplots*(Nfig-1))
          PlN = (iPl-1)*numSubplots + iT; % plotnumber
          cfgtopo.xlim = [time(ind_timewin_min+PlN-1) time(ind_timewin_min+PlN-1)];
          cfgtopo.highlightchannel   = list{PlN};
          if hastime
            cfgtopo.comment = strcat('time: ',num2str(time(ind_timewin_min+PlN-1)), ' s');
          elseif hasfreq
            cfgtopo.comment = strcat('freq: ',num2str(time(ind_timewin_min+PlN-1)), ' Hz');
          end
          cfgtopo.commentpos = 'title';
          subplot(cfg.subplotsize(1), cfg.subplotsize(2), iT);
          count = count+1;
          fprintf('making subplot %d from %d\n', count, Npl);
          ft_topoplotTFR(cfgtopo, stat);
        end
      end
    else
      cfgtopo.highlightchannel = list{1};
      cfgtopo.comment = strcat(compos, comneg);
      cfgtopo.commentpos = 'title';
      count = count+1;
      fprintf('making subplot %d from %d\n', count, Npl);
      ft_topoplotTFR(cfgtopo, stat);
    end
    if isequal(cfg.saveaspng, 'no')
      % nothing to do
    else
      % save figure
      filename = strcat(cfg.saveaspng, '_fig', num2str(iPl));
      print(gcf, '-dpng', filename);
    end
  end
end

ft_progress('close');

% return to previous warning settings
ft_warning(ws);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous stat
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sign = prob2hlsign(prob, hlsign)
if prob < 0.01
  sign = hlsign(1);
elseif prob < 0.05
  sign = hlsign(2);
elseif prob < 0.1
  sign = hlsign(3);
elseif prob < 0.2
  sign = hlsign(4);
elseif prob < 0.3
  sign = hlsign(5);
end
