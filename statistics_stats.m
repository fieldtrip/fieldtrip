function [stat, cfg] = statistics_stats(cfg, dat, design);

% This is a helper function that performs a massive univariate statistical
% test. This function is called by either TIMELOCKSTATISTICS, FREQSTATISTICS
% or SOURCSTATISTICS.
%
%  This function uses the Matlab statistics toolbox to perform various
%  statistical tests on timelock, frequency or source data. Supported
%  configuration options are
%   cfg.alpha     = 0.05
%   cfg.tail      = 0, -1 or 1
%   cfg.feedback  = 'no', 'text', 'textbar', 'gui'
%   cfg.method    = 'stats'
%   cfg.statistic = 'ttest'        test against a mean of zero
%                   'ttest2'	     compare the mean in two conditions
%                   'paired-ttest'
%                   'anova1'
%                   'kruskalwallis'
%
% See also TTEST, TTEST2, KRUSKALWALLIS

% Undocumented local options:
% cfg.avgovertime
% cfg.constantvalue

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: statistics_stats.m,v $
% Revision 1.8  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.7  2006/06/07 12:59:09  roboos
% fixed a bug in determining the size of the design matrix
%
% Revision 1.6  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.5  2006/04/10 12:13:15  roboos
% besides prob, also return mask and stat when possible
% return cfg as second output argument
%
% Revision 1.4  2006/03/30 12:24:34  roboos
% Implemented private/fixinside, which facilitates consistent
% handling of source/volume data. Improved documentation. Fixed some
% bugs related to inconsistent handling of ROIs (i.e. inside/outside)
%
% Revision 1.3  2005/12/08 16:58:20  ingnie
% fixed bug paired ttest, added log
%

fieldtripdefs

% test for the presence of the statistics toolbox
hasstats = (exist('ttest') & exist('ttest2'));
if ~hasstats
  error('this function requires the Matlab statistics toolbox');
end

% set the defaults that are common to all methods
if ~isfield(cfg, 'feedback'), cfg.feedback = 'textbar'; end

switch cfg.statistic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'ttest', 'ttest_samples_vs_const'}

  % set the defaults
  if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  if ~isfield(cfg, 'constantvalue'), cfg.constantvalue = 0; end
  if ~isfield(cfg, 'tail'), cfg.tail = 0; end

  if ~any(size(design)==1)
    error('design matrix should only contain one factor (i.e. one row)');
  end
  Ncond = length(unique(design));
  if Ncond>1
    error(sprintf('%s method is only supported for one condition at a time', cfg.statistic));
  end
  Nobs  = size(dat, 1);
  Nrepl = size(dat, 2); % over all conditions

  h = zeros(Nobs, 1);
  p = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);

  progress('init', cfg.feedback);
  for chan = 1:Nobs
    progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :)] = ttest(dat(chan, :), cfg.constantvalue, cfg.alpha, cfg.tail);
  end
  progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'ttest2', 'ttest_2samples_by_timepoint'}

  % set the defaults
  if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  if ~isfield(cfg, 'tail'), cfg.tail = 0; end

  if size(design,1)~=1
    error('design matrix should only contain one factor (i.e. one row)');
  end
  Ncond = length(unique(design));
  if Ncond~=2
    error(sprintf('%s method is only supported for two condition', cfg.statistic));
  end
  Nobs  = size(dat, 1);
  selA = find(design==design(1));
  selB = find(design~=design(1));
  Nrepl = [length(selA), length(selB)];

  h = zeros(Nobs, 1);
  p = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  progress('init', cfg.feedback);
  for chan = 1:Nobs
    progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :), stats] = ttest2(dat(chan, selA), dat(chan, selB), cfg.alpha, cfg.tail);
    s(chan) = stats.tstat;
  end
  progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'paired-ttest'}

  % set the defaults
  if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  if ~isfield(cfg, 'tail'), cfg.tail = 0; end

  if ~any(size(design)==1)
    error('design matrix should only contain one factor (i.e. one row)');
  end
  Ncond = length(unique(design));
  if Ncond~=2
    error(sprintf('%s method is only supported for two condition', cfg.statistic));
  end
  Nobs  = size(dat, 1);
  selA = find(design==design(1));
  selB = find(design~=design(1));
  Nrepl = [length(selA), length(selB)];
  if Nrepl(1)~=Nrepl(2)
    error('number of replications per condition should be the same');
  end

  h = zeros(Nobs, 1);
  p = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  progress('init', cfg.feedback);
  for chan = 1:Nobs
    progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :)] = ttest(dat(chan, selA)-dat(chan, selB), 0, cfg.alpha, cfg.tail);
  end
  progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'anova1'}

  if ~any(size(design)==1)
    error('design matrix should only contain one factor (i.e. one row)');
  end
  Ncond = length(unique(design));
  Nobs  = size(dat, 1);
  Nrepl = size(dat, 2); % over all conditions

  h = zeros(Nobs, 1);
  p = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);
  fprintf('number of levels %d\n', Ncond);

  progress('init', cfg.feedback);
  for chan = 1:Nobs
    progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    p(chan) = anova1(dat(chan, :), design(:), 'off');
  end
  progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'kruskalwallis'}

  if ~any(size(design)==1)
    error('design matrix should only contain one factor (i.e. one row)');
  end
  Ncond = length(unique(design));
  Nobs  = size(dat, 1);
  Nrepl = size(dat, 2); % over all conditions

  h = zeros(Nobs, 1);
  p = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);
  fprintf('number of levels %d\n', Ncond);

  progress('init', cfg.feedback);
  for chan = 1:Nobs
    progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    p(chan) = kruskalwallis(dat(chan, :), design(:), 'off');
  end
  progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case {'anovan'}
% 
%   Nfact = size(design,1);
%   Nobs  = size(dat, 1);
%   Nrepl = size(dat, 2); % over all conditions
% 
%   h = zeros(Nobs, 1);
%   p = zeros(Nobs, 1);
%   ci = zeros(Nobs, 2);
%   fprintf('number of observations %d\n', Nobs);
%   fprintf('number of replications %d\n', Nrepl);
%   fprintf('number of factors %d\n', Nfact);
% 
%   % reformat the design matrix into the grouping variable cell-array
%   for i=1:Nfact
%     group{i} = design(i,:);
%   end
% 
%   progress('init', cfg.feedback);
%   for chan = 1:Nobs
%     progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
%     % FIXME, the probability is returned for each factor separately
%     p = anovan(dat(chan, :), group, 'display', 'off');
%   end
%   progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ttest_window_avg_vs_const'
  % this used to be a feature of the timelockanaolysis as it was
  % originally implemented by Jens Schwartzbach, but it has been
  % superseded by the use of prepare_timefreq_data for data selection
  error(sprintf('%s is not supported any more, use cfg.avgovertime=''yes'' instead', cfg.statistic));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
  error(sprintf('Statistical method ''%s'' is not implemented', cfg.statistic));
end

% assign the output variable
stat = [];
try, stat.mask = h; end
try, stat.prob = p; end
try, stat.stat = s; end

