function [stat, cfg] = ft_statistics_stats(cfg, dat, design)

% FT_STATISTICS_STATS performs a massive univariate statistical test using the
% MATLAB statistics toolbox. This function should not be called directly,
% instead you should call the function that is associated with the type of data
% on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% Where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS
% or FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively and with
% cfg.method = 'montecarlo'
%
%  This function uses the MATLAB statistics toolbox to perform various
%  statistical tests on timelock, frequency or source data. Supported
%  configuration options are
%   cfg.alpha     = number, critical value for rejecting the null-hypothesis (default = 0.05)
%   cfg.tail      = number, -1, 1 or 0 (default = 0)
%   cfg.feedback  = string, 'gui', 'text', 'textbar' or 'no' (default = 'textbar')
%   cfg.method    = 'stats'
%   cfg.statistic = 'ttest'          test against a mean of zero
%                   'ttest2'         compare the mean in two conditions
%                   'paired-ttest'
%                   'anova1'
%                   'kruskalwallis'
%                   'signtest'
%                   'signrank'
%
% See also TTEST, TTEST2, KRUSKALWALLIS, SIGNTEST, SIGNRANK

% Undocumented local options:
% cfg.avgovertime
% cfg.constantvalue

% Copyright (C) 2005, Robert Oostenveld
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

% test for the presence of the statistics toolbox
ft_hastoolbox('stats', 1);

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
  s = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :), stats] = ttest_wrapper(dat(chan, :), cfg.constantvalue, cfg.alpha, cfg.tail);
    s(chan) = stats.tstat;
  end
  ft_progress('close');

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
  s = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :), stats] = ttest2_wrapper(dat(chan, selA), dat(chan, selB), cfg.alpha, cfg.tail);
    s(chan) = stats.tstat;
  end
  ft_progress('close');

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
  s = zeros(Nobs, 1);
  ci = zeros(Nobs, 2);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [h(chan), p(chan), ci(chan, :), stats] = ttest_wrapper(dat(chan, selA)-dat(chan, selB), 0, cfg.alpha, cfg.tail);
    s(chan) = stats.tstat;
  end
  ft_progress('close');

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
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);
  fprintf('number of levels %d\n', Ncond);

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    p(chan) = anova1(dat(chan, :), design(:), 'off');
  end
  ft_progress('close');

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
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d\n', Nrepl);
  fprintf('number of levels %d\n', Ncond);

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    p(chan) = kruskalwallis(dat(chan, :), design(:), 'off');
  end
  ft_progress('close');

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
%   ft_progress('init', cfg.feedback);
%   for chan = 1:Nobs
%     ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
%     % FIXME, the probability is returned for each factor separately
%     p = anovan(dat(chan, :), group, 'display', 'off');
%   end
%   ft_progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ttest_window_avg_vs_const'
  % this used to be a feature of the timelockanalysis as it was
  % originally implemented by Jens Schwartzbach, but it has been
  % superseded by the use of ft_selectdata for data selection
  error(sprintf('%s is not supported any more, use cfg.avgovertime=''yes'' instead', cfg.statistic));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'signtest'}
  % set the defaults
  if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  if ~isfield(cfg, 'tail'), cfg.tail = 0; end
  
  switch cfg.tail
    case 0
      cfg.tail = 'both';
    case -1
      cfg.tail = 'left';
    case 1
      cfg.tail = 'right';
  end;
  
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
  s = zeros(Nobs, 1);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [p(chan), h(chan), stats] = signtest(dat(chan, selA), dat(chan, selB),'alpha', cfg.alpha,'tail', cfg.tail);
    s(chan) = stats.sign;
  end
  ft_progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'signrank'}
  % set the defaults
  if ~isfield(cfg, 'alpha'), cfg.alpha = 0.05; end
  if ~isfield(cfg, 'tail'), cfg.tail = 0; end
  
  switch cfg.tail
    case 0
      cfg.tail = 'both';
    case -1
      cfg.tail = 'left';
    case 1
      cfg.tail = 'right';
  end;
  
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
  s = zeros(Nobs, 1);
  fprintf('number of observations %d\n', Nobs);
  fprintf('number of replications %d and %d\n', Nrepl(1), Nrepl(2));

  ft_progress('init', cfg.feedback);
  for chan = 1:Nobs
    ft_progress(chan/Nobs, 'Processing observation %d/%d\n', chan, Nobs);
    [p(chan), h(chan), stats] = signrank(dat(chan, selA), dat(chan, selB),'alpha', cfg.alpha,'tail', cfg.tail);
    s(chan) = stats.signedrank;
  end
  ft_progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
  error(sprintf('Statistical method ''%s'' is not implemented', cfg.statistic));
end

% assign the output variable
stat = [];
try, stat.mask = h; end
try, stat.prob = p; end
try, stat.stat = s; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions for ttest and ttest2
% - old Matlab, with syntax:                   ttest(x,y,alpha,tail,dim)
% - new Matlab and GNU Octave, with syntax:    ttest(x,y,'alpha',alpha,...)
function [h,p,ci,stats]=ttest_wrapper(x,y,alpha,tail)
    [h,p,ci,stats]=general_ttestX_wrapper(@ttest,x,y,alpha,tail);

function [h,p,ci,stats]=ttest2_wrapper(x,y,alpha,tail)
    [h,p,ci,stats]=general_ttestX_wrapper(@ttest2,x,y,alpha,tail);

function [h,p,ci,stats]=general_ttestX_wrapper(ttest_func,x,y,alpha,tail)
    if nargin(ttest_func)>0
        % old Matlab
        [h,p,ci,stats]=ttest_func(x,y,alpha,tail);

    else
        % GNU Octave and new Matlab
        [h,p,ci,stats]=ttest_func(x,y,'alpha',alpha,'tail',tail);
    end

