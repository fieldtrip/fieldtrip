function [cfg] = ft_clusterplot(cfg, stat)

% FT_CLUSTERPLOT plots a series of topoplots with found clusters highlighted.
% stat is 2D or 1D data from FT_TIMELOCKSTATISTICS or FT_FREQSTATISTICS with 'cluster'
% as cfg.correctmc. 2D: stat from FT_TIMELOCKSTATISTICS not averaged over
% time, or stat from FT_FREQSTATISTICS averaged over frequency not averaged over
% time. 1D: averaged over time as well.
%
% Use as
%   ft_clusterplot(cfg, stat)
%
% Where the configuration options can be
%   cfg.alpha                     = number, highest cluster p-value to be plotted
%                                   max 0.3 (default = 0.05)
%   cfg.highlightseries           = 1x5 cell-array, highlight option series ('on','labels','numbers')
%                                   default {'on','on','on','on','on'} for p < [0.01 0.05 0.1 0.2 0.3]
%   cfg.highlightsymbolseries     = 1x5 vector, highlight marker symbol series
%                                   default ['*','x','+','o','.'] for p < [0.01 0.05 0.1 0.2 0.3]
%   cfg.highlightsizeseries       = 1x5 vector, highlight marker size series
%                                   default [6 6 6 6 6] for p < [0.01 0.05 0.1 0.2 0.3]
%   cfg.highlightcolorpos         = color of highlight marker for positive clusters
%                                   default = [0 0 0]
%   cfg.highlightcolorneg         = color of highlight marker for negative clusters
%                                   default = [0 0 0]
%   cfg.saveaspng                 = string, filename of the output figures (default = 'no')
%
% It is also possible to specify other cfg options that apply to FT_TOPOPLOTTFR.
% You CANNOT specify cfg.xlim, any of the FT_TOPOPLOTTFR highlight
% options, cfg.comment and cfg.commentpos.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also:
%   FT_TOPOPLOTTFR, FT_TOPOPLOTER, FT_SINGLEPLOTER

% Copyright (C) 2007, Ingrid Nieuwenhuis, F.C. Donders Centre
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar stat

% check if the given data is appropriate
if isfield(stat,'freq') && length(stat.freq) > 1
  error('stat contains multiple frequencies which is not allowed because it should be averaged over frequencies')
end

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
if ~isfield(cfg,'alpha'),                  cfg.alpha = 0.05;                                    end;
if ~isfield(cfg,'highlightseries'),        cfg.highlightseries = {'on','on','on','on','on'};    end;
if ~isfield(cfg,'highlightsymbolseries'),  cfg.highlightsymbolseries = ['*','x','+','o','.'];   end;
if ~isfield(cfg,'highlightsizeseries'),    cfg.highlightsizeseries = [6 6 6 6 6];               end;
if ~isfield(cfg,'hllinewidthseries'),      cfg.hllinewidthseries = [1 1 1 1 1];                 end;
if ~isfield(cfg,'highlightcolorpos'),      cfg.highlightcolorpos = [0 0 0];                     end;
if ~isfield(cfg,'highlightcolorneg'),      cfg.highlightcolorneg = [0 0 0];                     end;
if ~isfield(cfg,'parameter'),              cfg.parameter = 'stat';                              end;
if ~isfield(cfg,'saveaspng'),              cfg.saveaspng = 'no';                                end;

% error if cfg.highlightseries is not a cell, for possible confusion with cfg-options
if ~iscell(cfg.highlightseries)
  error('cfg.highlightseries should be a cell-array of strings')
end

% set additional options for topoplotting
if isfield(cfg, 'marker'),                cfgtopo.marker         = cfg.marker ;         end
if ~isfield(cfg,'marker'),                cfgtopo.marker         = 'off';               end
if isfield(cfg, 'markersymbol'),          cfgtopo.markersymbol   = cfg.markersymbol;    end
if isfield(cfg, 'markercolor'),           cfgtopo.markercolor    = cfg.markercolor;     end
if isfield(cfg, 'markersize'),            cfgtopo.markersize     = cfg.markersize;      end
if isfield(cfg, 'markerfontsize'),        cfgtopo.markerfontsize = cfg.markerfontsize;  end
if isfield(cfg, 'style'),                 cfgtopo.style          = cfg.style ;          end
if isfield(cfg, 'gridscale'),             cfgtopo.gridscale      = cfg.gridscale;       end
if isfield(cfg, 'interplimits'),          cfgtopo.interplimits   = cfg.interplimits;    end
if isfield(cfg, 'interpolation'),         cfgtopo.interpolation  = cfg.interpolation;   end
if isfield(cfg, 'contournum'),            cfgtopo.contournum     = cfg.contournum;      end
if isfield(cfg, 'colorbar'),              cfgtopo.colorbar       = cfg.colorbar;        end
if isfield(cfg, 'shading'),               cfgtopo.shading        = cfg.shading';        end
if isfield(cfg, 'zlim'),                  cfgtopo.zlim           = cfg.zlim;            end

cfgtopo.parameter = cfg.parameter;

% prepare the layout, this only has to be done once
cfgtopo.layout = ft_prepare_layout(cfg, stat);

% detect 2D or 1D
is2D = isfield(stat,'time');

% add .time field to 1D data, topoplotER wants it
if ~is2D
  stat.time = 0; %doesn't matter what it is, so just choose 0
end;

% find significant clusters
sigpos = [];
signeg = [];
haspos = isfield(stat,'posclusters');
hasneg = isfield(stat,'negclusters');

if haspos == 0 && hasneg == 0
  fprintf('%s\n','no significant clusters in data; nothing to plot')
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
  sigpos = find(sigpos == 1);
  signeg = find(signeg == 1);
  Nsigpos = length(sigpos);
  Nsigneg = length(signeg);
  Nsigall = Nsigpos + Nsigneg;
  
  if Nsigall == 0
    error('no clusters present with a p-value lower than the specified alpha, nothing to plot')
  end
  
  % make clusterslabel matrix per significant cluster
  if haspos
    posCLM = squeeze(stat.posclusterslabelmat);
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
    negCLM = squeeze(stat.negclusterslabelmat);
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
  
  fprintf('%s%i%s%g%s\n','There are ',Nsigall,' clusters smaller than alpha (',cfg.alpha,')')
  
  if is2D
    % define time window per cluster
    for iPos = 1:length(sigpos)
      possum_perclus = sum(sigposCLM(:,:,iPos),1); %sum over Chans for each timepoint
      ind_min = min(find(possum_perclus~=0));
      ind_max = max(find(possum_perclus~=0));
      time_perclus = [stat.time(ind_min) stat.time(ind_max)];
      fprintf('%s%s%s%s%s%s%s%s%s%s%s\n','Positive cluster: ',num2str(sigpos(iPos)),', pvalue: ',num2str(probpos(iPos)),' (',hlsignpos(iPos),')',', t = ',num2str(time_perclus(1)),' to ',num2str(time_perclus(2)))
    end
    for iNeg = 1:length(signeg)
      negsum_perclus = sum(signegCLM(:,:,iNeg),1);
      ind_min = min(find(negsum_perclus~=0));
      ind_max = max(find(negsum_perclus~=0));
      time_perclus = [stat.time(ind_min) stat.time(ind_max)];
      fprintf('%s%s%s%s%s%s%s%s%s%s%s\n','Negative cluster: ',num2str(signeg(iNeg)),', pvalue: ',num2str(probneg(iNeg)),' (',hlsignneg(iNeg),')',', t = ',num2str(time_perclus(1)),' to ',num2str(time_perclus(2)))
    end
    
    % define timewindow containing all significant clusters
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
    
    ind_timewin_min = min(find(allsum~=0));
    ind_timewin_max = max(find(allsum~=0));
    timewin = stat.time(ind_timewin_min:ind_timewin_max);
    
  else
    for iPos = 1:length(sigpos)
      fprintf('%s%s%s%s%s%s%s\n','Positive cluster: ',num2str(sigpos(iPos)),', pvalue: ',num2str(probpos(iPos)),' (',hlsignpos(iPos),')')
    end
    for iNeg = 1:length(signeg)
      fprintf('%s%s%s%s%s%s%s\n','Negative cluster: ',num2str(signeg(iNeg)),', pvalue: ',num2str(probneg(iNeg)),' (',hlsignneg(iNeg),')')
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
    elseif stat.posclusters(sigpos(iPos)).prob < 0.05
      cfgtopo.highlight{iPos}         = cfg.highlightseries{2};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(2);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(2);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.1
      cfgtopo.highlight{iPos}         = cfg.highlightseries{3};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(3);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(3);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.2
      cfgtopo.highlight{iPos}         = cfg.highlightseries{4};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(4);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(4);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.3
      cfgtopo.highlight{iPos}         = cfg.highlightseries{5};
      cfgtopo.highlightsymbol{iPos}   = cfg.highlightsymbolseries(5);
      cfgtopo.highlightsize{iPos}     = cfg.highlightsizeseries(5);
    end
    cfgtopo.highlightcolor{iPos}        = cfg.highlightcolorpos;
    compos = strcat(compos,cfgtopo.highlightsymbol{iPos}, 'p=',num2str(probpos(iPos)),' '); % make comment, only used for 1D data
  end
  
  for iNeg = 1:length(signeg)
    if stat.negclusters(signeg(iNeg)).prob < 0.01
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{1};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(1);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(1);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.05
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{2};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(2);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(2);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.1
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{3};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(3);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(3);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.2
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{4};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(4);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(4);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.3
      cfgtopo.highlight{length(sigpos)+iNeg}         = cfg.highlightseries{5};
      cfgtopo.highlightsymbol{length(sigpos)+iNeg}   = cfg.highlightsymbolseries(5);
      cfgtopo.highlightsize{length(sigpos)+iNeg}     = cfg.highlightsizeseries(5);
    end
    cfgtopo.highlightcolor{length(sigpos)+iNeg}        = cfg.highlightcolorneg;
    comneg = strcat(comneg,cfgtopo.highlightsymbol{length(sigpos)+iNeg}, 'p=',num2str(probneg(iNeg)),' '); % make comment, only used for 1D data
  end
  
  if is2D
    Npl = length(timewin);
  else
    Npl = 1;
  end
  Nfig = ceil(Npl/15);
  
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
  
  
  % make plots
  for iPl = 1:Nfig
    figure;
    if is2D
      if iPl < Nfig
        for iT = 1:15
          PlN = (iPl-1)*15 + iT; %plotnumber
          cfgtopo.xlim = [stat.time(ind_timewin_min+PlN-1) stat.time(ind_timewin_min+PlN-1)];
          cfgtopo.highlightchannel = list{PlN};
          cfgtopo.comment = strcat('time: ',num2str(stat.time(ind_timewin_min+PlN-1)), ' s');
          cfgtopo.commentpos = 'title';
          subplot(3,5,iT);
          ft_topoplotTFR(cfgtopo, stat);
        end
      elseif iPl == Nfig
        for iT = 1:Npl-(15*(Nfig-1))
          PlN = (iPl-1)*15 + iT; %plotnumber
          cfgtopo.xlim = [stat.time(ind_timewin_min+PlN-1) stat.time(ind_timewin_min+PlN-1)];
          cfgtopo.highlightchannel   = list{PlN};
          cfgtopo.comment = strcat('time: ',num2str(stat.time(ind_timewin_min+PlN-1)), ' s');
          cfgtopo.commentpos = 'title';
          subplot(3,5,iT);
          ft_topoplotTFR(cfgtopo, stat);
        end
      end
    else
      cfgtopo.highlightchannel = list{1};
      cfgtopo.xparam = 'time';
      cfgtopo.yparam = '';
      cfgtopo.comment = strcat(compos,comneg);
      cfgtopo.commentpos = 'title';
      ft_topoplotTFR(cfgtopo, stat);
    end
    % save figure
    if isequal(cfg.saveaspng,'no');
    else
      filename = strcat(cfg.saveaspng, '_fig', num2str(iPl));
      print(gcf,'-dpng',filename);
    end
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous stat

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

