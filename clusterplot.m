function clusterplot(cfg, stat)

% CLUSTERPLOT plots a series of topoplots with found clusters highlighted.
% stat is 2D or 1D data from TIMELOCKSTATISTICS or FREQSTATISTICS with 'cluster'
% as cfg.correctmc. 2D: stat from timelockstatistics not averaged over
% time, or stat from freqstatistics averaged over frequency not averaged over
% time. 1D: averaged over time as well.
%
% use as: clusterplot(cfg,stat)
%
% configuration options
% cfg.alpha              = number, highest cluster p-value to be plotted
%                          max 0.3 (default = 0.05)
% cfg.hlmarkerseries     = 1x5 vector, highlight marker symbol series
%                          default ['*','x','+','o','.'] for p < [0.01 0.05 0.1 0.2 0.3]
% cfg.hlmarkersizeseries = 1x5 vector, highlight marker size series
%                          default [6 6 6 6 6] for p < [0.01 0.05 0.1 0.2 0.3]
% cfg.hllinewidthseries  = 1x5 vector, highlight marker linewidth series
%                          default [1 1 1 1 1] for p < [0.01 0.05 0.1 0.2 0.3]
% cfg.hlcolorpos         = color of highlight marker for positive clusters
%                          default = [0 0 0]
% cfg.hlcolorneg         = color of highlight marker for negative clusters
%                          default = [0 0 0]
% cfg.saveaspng          = string, path where figure has to be saved to (default = 'no')
%                          When multiple figures figure gets extension with fignum
%
% It is also possible to specify other cfg options that apply to TOPOPLOTER
% or TOPOPLOT. You CANNOT specify cfg.xlim, any of the TOPOPLOT highlight
% options, cfg.comment and cfg.commentpos.
%
% See also:
%   topoplot, singleplotER

% Copyright (C) 2007, Ingrid Nieuwenhuis, F.C. Donders Centre
%
% $Log: clusterplot.m,v $
% Revision 1.8  2009/05/22 14:54:58  ingnie
% added check of data, should be averged over frequencies
%
% Revision 1.7  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.6  2008/09/22 15:18:53  roboos
% only prepare the layout once
%
% Revision 1.5  2008/06/12 12:22:36  ingnie
% Now also works for 1D data (= avaraged over time), added some comma's in default setting for readability
%
% Revision 1.4  2007/11/07 16:28:08  ingnie
% added option saveaspng
%
% Revision 1.3  2007/11/07 12:52:46  ingnie
% add cvs log
%

fieldtripdefs

% check if given data is appropriate
if isfield(stat,'freq') && length(stat.freq) > 1
  error('stat contains multiple frequencies which is not allowed because it should be averaged over frequencies')
end

% set the defaults
if ~isfield(cfg,'alpha'),                  cfg.alpha = 0.05;                             end;
if ~isfield(cfg,'hlmarkerseries'),         cfg.hlmarkerseries = ['*','x','+','o','.'];   end;
if ~isfield(cfg,'hlmarkersizeseries'),     cfg.hlmarkersizeseries = [6 6 6 6 6];         end;
if ~isfield(cfg,'hllinewidthseries'),      cfg.hllinewidthseries = [1 1 1 1 1];          end;
if ~isfield(cfg,'hlcolorpos'),             cfg.hlcolorpos = [0 0 0];                     end;
if ~isfield(cfg,'hlcolorneg'),             cfg.hlcolorneg = [0 0 0];                     end;
if ~isfield(cfg,'zparam'),                 cfg.zparam = 'stat';                          end;
if ~isfield(cfg,'saveaspng'),              cfg.saveaspng = 'no';                         end;

% prepare the layout, this only has to be done once
cfg.layout = prepare_layout(cfg, stat);

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
  for iPos = 1:length(stat.posclusters)
    sigpos(iPos) = stat.posclusters(iPos).prob < cfg.alpha;
  end
  for iNeg = 1:length(stat.negclusters)
    signeg(iNeg) = stat.negclusters(iNeg).prob < cfg.alpha;
  end
  sigpos = find(sigpos == 1);
  signeg = find(signeg == 1);
  Nsigpos = length(sigpos);
  Nsigneg = length(signeg);
  Nsigall = Nsigpos + Nsigneg;

  % make clusterslabel matrix per significant cluster
  posCLM = squeeze(stat.posclusterslabelmat);
  sigposCLM = zeros(size(posCLM));
  probpos = [];
  for iPos = 1:length(sigpos)
    sigposCLM(:,:,iPos) = (posCLM == sigpos(iPos));
    probpos(iPos) = stat.posclusters(iPos).prob;
    hlsignpos(iPos) = prob2hlsign(probpos(iPos), cfg.hlmarkerseries);
  end

  negCLM = squeeze(stat.negclusterslabelmat);
  signegCLM = zeros(size(negCLM));
  probneg = [];
  for iNeg = 1:length(signeg)
    signegCLM(:,:,iNeg) = (negCLM == signeg(iNeg));
    probneg(iNeg) = stat.negclusters(iNeg).prob;
    hlsignneg(iNeg) = prob2hlsign(probneg(iNeg), cfg.hlmarkerseries);
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
    allsum = possum + negsum;

    ind_timewin_min = min(find(allsum~=0));
    ind_timewin_max = max(find(allsum~=0));

    timestep = stat.time(2) - stat.time(1);
    timewin = [stat.time(ind_timewin_min): timestep :stat.time(ind_timewin_max)];
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
      cfg.hlmarker{iPos}     = cfg.hlmarkerseries(1);
      cfg.hlmarkersize{iPos} = cfg.hlmarkersizeseries(1);
      cfg.hllinewidth{iPos}  = cfg.hllinewidthseries(1);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.05
      cfg.hlmarker{iPos}     = cfg.hlmarkerseries(2);
      cfg.hlmarkersize{iPos} = cfg.hlmarkersizeseries(2);
      cfg.hllinewidth{iPos}  = cfg.hllinewidthseries(2);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.1
      cfg.hlmarker{iPos}     = cfg.hlmarkerseries(3);
      cfg.hlmarkersize{iPos} = cfg.hlmarkersizeseries(3);
      cfg.hllinewidth{iPos}  = cfg.hllinewidthseries(3);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.2
      cfg.hlmarker{iPos}     = cfg.hlmarkerseries(4);
      cfg.hlmarkersize{iPos} = cfg.hlmarkersizeseries(4);
      cfg.hllinewidth{iPos}  = cfg.hllinewidthseries(4);
    elseif stat.posclusters(sigpos(iPos)).prob < 0.3
      cfg.hlmarker{iPos}     = cfg.hlmarkerseries(5);
      cfg.hlmarkersize{iPos} = cfg.hlmarkersizeseries(5);
      cfg.hllinewidth{iPos}  = cfg.hllinewidthseries(5);
    end
    cfg.hlcolor{iPos}        = cfg.hlcolorpos;
    compos = strcat(compos,cfg.hlmarker{iPos}, 'p=',num2str(probpos(iPos)),' '); % make comment, only used for 1D data
  end

  for iNeg = 1:length(signeg)
    if stat.negclusters(signeg(iNeg)).prob < 0.01
      cfg.hlmarker{length(sigpos)+iNeg}     = cfg.hlmarkerseries(1);
      cfg.hlmarkersize{length(sigpos)+iNeg} = cfg.hlmarkersizeseries(1);
      cfg.hllinewidth{length(sigpos)+iNeg}  = cfg.hllinewidthseries(1);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.05
      cfg.hlmarker{length(sigpos)+iNeg}     = cfg.hlmarkerseries(2);
      cfg.hlmarkersize{length(sigpos)+iNeg} = cfg.hlmarkersizeseries(2);
      cfg.hllinewidth{length(sigpos)+iNeg}  = cfg.hllinewidthseries(2);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.1
      cfg.hlmarker{length(sigpos)+iNeg}     = cfg.hlmarkerseries(3);
      cfg.hlmarkersize{length(sigpos)+iNeg} = cfg.hlmarkersizeseries(3);
      cfg.hllinewidth{length(sigpos)+iNeg}  = cfg.hllinewidthseries(3);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.2
      cfg.hlmarker{length(sigpos)+iNeg}     = cfg.hlmarkerseries(4);
      cfg.hlmarkersize{length(sigpos)+iNeg} = cfg.hlmarkersizeseries(4);
      cfg.hllinewidth{length(sigpos)+iNeg}  = cfg.hllinewidthseries(4);
    elseif stat.negclusters(signeg(iNeg)).prob < 0.3
      cfg.hlmarker{length(sigpos)+iNeg}     = cfg.hlmarkerseries(5);
      cfg.hlmarkersize{length(sigpos)+iNeg} = cfg.hlmarkersizeseries(5);
      cfg.hllinewidth{length(sigpos)+iNeg}  = cfg.hllinewidthseries(5);
    end
    cfg.hlcolor{length(sigpos)+iNeg}        = cfg.hlcolorneg;
    comneg = strcat(comneg,cfg.hlmarker{length(sigpos)+iNeg}, 'p=',num2str(probneg(iNeg)),' '); % make comment, only used for 1D data
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
          cfg.xlim = [stat.time(ind_timewin_min+PlN-1) stat.time(ind_timewin_min+PlN-1)];
          cfg.highlight = list{PlN};
          cfg.comment = strcat('time: ',num2str(stat.time(ind_timewin_min+PlN-1)), ' s');
          cfg.commentpos = 'title';
          subplot(3,5,iT);
          topoplotER(cfg, stat);
        end
      elseif iPl == Nfig
        for iT = 1:Npl-(15*(Nfig-1))
          PlN = (iPl-1)*15 + iT; %plotnumber
          cfg.xlim = [stat.time(ind_timewin_min+PlN-1) stat.time(ind_timewin_min+PlN-1)];
          cfg.highlight   = list{PlN};
          cfg.comment = strcat('time: ',num2str(stat.time(ind_timewin_min+PlN-1)), ' s');
          cfg.commentpos = 'title';
          subplot(3,5,iT);
          topoplotER(cfg, stat);
        end
      end
    else
      cfg.highlight = list{1};
      cfg.xparam = 'time';
      cfg.yparam = '';
      cfg.comment = strcat(compos,comneg);
      cfg.commentpos = 'title';
      topoplotER(cfg, stat);
    end
    % save figure
    if isequal(cfg.saveaspng,'no');
    else
      filename = strcat(cfg.saveaspng, '_fig', num2str(iPl));
      print(gcf,'-dpng',filename);
    end
  end
end

%% subfunctions %%
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
