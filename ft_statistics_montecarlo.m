function [stat, cfg] = ft_statistics_montecarlo(cfg, dat, design, varargin)

% FT_STATISTICS_MONTECARLO performs a nonparametric statistical test by calculating
% Monte-Carlo estimates of the significance probabilities and/or critical values from
% the permutation distribution. It is not recommended to call this function directly,
% but instead you should call the function that is associated with the type of data on 
% which you want to perform the test. This is because important data bookkeeping 
% operations on the data are performed in the higher-level functions, which are
% assumed to have been handled correctly for the input arguments into this function. 
% Also, notably, a prespecified randomseed in the cfg is handled in the higher
% level function, not here.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS or
% FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively 
% and with cfg.method = 'montecarlo'. 
%
% The configuration options that can be specified are:
%   cfg.numrandomization = number of randomizations, can be 'all'
%   cfg.correctm         = string, apply multiple-comparison correction, 'no', 'max', cluster', 'tfce', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
%   cfg.alpha            = number, critical value for rejecting the null-hypothesis per tail (default = 0.05)
%   cfg.tail             = number, -1, 1 or 0 (default = 0)
%   cfg.correcttail      = string, correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s)
%   cfg.wvar             = number or list with indices, within-cell variable(s)
%   cfg.cvar             = number or list with indices, control variable(s)
%   cfg.feedback         = string, 'gui', 'text', 'textbar' or 'no' (default = 'text')
%   cfg.randomseed       = string, 'yes', 'no' or a number (default = 'yes'), this option is not used in this 
%                            function directly, but is handled by ft_<something>statistics. If you want to control
%                            the random number generator while calling ft_statistics_montecarlo directly, it should 
%                            be specified on the command line (or in the caller script) just prior to calling this
%                            function.
%
% If you use a cluster-based statistic, you can specify the following options that
% determine how the single-sample or single-voxel statistics will be thresholded and
% combined into one statistical value per cluster.
%   cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%                          the option 'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity; 
%                          see Hayasaka & Nichols (2004) NeuroImage for details
%   cfg.clusterthreshold = method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%   cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
%   cfg.clustercritval   = for parametric thresholding (default is determined by the statfun)
%   cfg.clustertail      = -1, 1 or 0 (default = 0)
%
% To include the channel dimension for clustering of channel level data, you should specify
%   cfg.neighbours       = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
% If you specify an empty neighbourhood structure, clustering will only be done
% over frequency and/or time and not over neighbouring channels.
%
% The statistic that is computed for each sample in each random reshuffling
% of the data is specified as
%   cfg.statistic       = 'indepsamplesT'           independent samples T-statistic,
%                         'indepsamplesF'           independent samples F-statistic,
%                         'indepsamplesregrT'       independent samples regression coefficient T-statistic,
%                         'indepsamplesZcoh'        independent samples Z-statistic for coherence,
%                         'depsamplesT'             dependent samples T-statistic,
%                         'depsamplesFmultivariate' dependent samples F-statistic MANOVA,
%                         'depsamplesregrT'         dependent samples regression coefficient T-statistic,
%                         'actvsblT'                activation versus baseline T-statistic.
% or you can specify your own low-level statistical function.
%
% You can also use a custom statistic of your choice that is sensitive to the
% expected effect in the data. You can implement the statistic in a "statfun" that
% will be called for each randomization. The requirements on a custom statistical
% function is that the function is called ft_statfun_xxx, and that the function returns
% a structure with a "stat" field containing the single sample statistical values.
% Have a look at the functions in the fieldtrip/statfun directory (e.g. 
% FT_STATFUN_INDEPSAMPLEST) for the correct format of the input and output.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS,
% FT_STATISTICS_ANALYTIC, FT_STATISTICS_STATS, FT_STATISTICS_MVPA,
% FT_STATISTICS_CROSSVALIDATE

% Undocumented local options:
%   cfg.resampling       permutation, bootstrap
%   cfg.computecritval   yes|no, for the statfun
%   cfg.computestat      yes|no, for the statfun
%   cfg.computeprob      yes|no, for the statfun
%   cfg.voxelstatistic   deprecated
%   cfg.voxelthreshold   deprecated
%   cfg.precondition     before|after|[], for the statfun

% Copyright (C) 2005-2015, Robert Oostenveld
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

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

% check whether the function has been called from ft_timelockstatistics, ft_freqstatistics, or ft_sourcestatistics
st = dbstack;
m  = mfilename;
if isscalar(st)
  ft_warning('It seems that %s has been called directly from the command line. This is not recommended, unless you know what you are doing', m);
elseif numel(st)>1 && ~ismember(st(2).name, {'ft_freqstatistics' 'ft_timelockstatistics' 'ft_sourcestatistics'})
  ft_warning('It seems that %s has not been called from one of the FT_XXXSTATISTICS functions. This is not recommended, unless you know what you are doing', m);
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'factor',           'ivar'});
cfg = ft_checkconfig(cfg, 'renamed',     {'unitfactor',       'uvar'});
cfg = ft_checkconfig(cfg, 'renamed',     {'repeatedmeasures', 'uvar'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'clusterthreshold', 'nonparametric', 'nonparametric_individual'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'yes', 'max'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'none', 'no'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'bonferoni', 'bonferroni'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'holms', 'holm'});
cfg = ft_checkconfig(cfg, 'required',    {'statistic'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'ztransform', 'removemarginalmeans', 'randomfactor', 'voxelthreshold', 'voxelstatistic'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'statfun', 'depsamplesF', 'ft_statfun_depsamplesFmultivariate'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'statfun', 'ft_statfun_depsamplesF', 'ft_statfun_depsamplesFmultivariate'});

% set the defaults for the main function
cfg.alpha        = ft_getopt(cfg, 'alpha',      0.05);
cfg.tail         = ft_getopt(cfg, 'tail',       0);
cfg.correctm     = ft_getopt(cfg, 'correctm',   'no');
cfg.resampling   = ft_getopt(cfg, 'resampling', 'permutation');
cfg.feedback     = ft_getopt(cfg, 'feedback',   'text');
cfg.ivar         = ft_getopt(cfg, 'ivar',       'all');
cfg.uvar         = ft_getopt(cfg, 'uvar',       []);
cfg.cvar         = ft_getopt(cfg, 'cvar',       []);
cfg.wvar         = ft_getopt(cfg, 'wvar',       []);
cfg.correcttail  = ft_getopt(cfg, 'correcttail',  'no');
cfg.precondition = ft_getopt(cfg, 'precondition', []);

% explicit check for option 'yes' in cfg.correctail.
if strcmp(cfg.correcttail, 'yes')
  ft_error('cfg.correcttail = ''yes'' is not allowed, use either ''prob'', ''alpha'' or ''no''')
end

if strcmp(cfg.correctm, 'tfce')
  % TODO this could require some better defaults
  cfg.connectivity = ft_getopt(cfg, 'connectivity', []);
  cfg.tfce_h0      = ft_getopt(cfg, 'tfce_h0', 0);
  cfg.tfce_H       = ft_getopt(cfg, 'tfce_H',  2);
  cfg.tfce_E       = ft_getopt(cfg, 'tfce_E',  0.5);
  cfg.tfce_nsteps  = ft_getopt(cfg, 'tfce_nsteps', 100);
else
  % these options only apply to tfce, to ensure appropriate configs they are forbidden when _not_ clustering
  cfg = ft_checkconfig(cfg, 'unused', {'tfce_h0', 'tfce_H', 'tfce_E', 'tfce_nsteps'});
end

if strcmp(cfg.correctm, 'cluster')
  % set the defaults for clustering
  cfg.clusterstatistic = ft_getopt(cfg, 'clusterstatistic', 'maxsum');
  cfg.clusterthreshold = ft_getopt(cfg, 'clusterthreshold', 'parametric');
  cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',     0.05);
  cfg.clustercritval   = ft_getopt(cfg, 'clustercritval',   []);
  cfg.clustertail      = ft_getopt(cfg, 'clustertail',      cfg.tail);
  cfg.connectivity     = ft_getopt(cfg, 'connectivity',     []); % the default is dealt with below
else
  % these options only apply to clustering, to ensure appropriate configs they are forbidden when _not_ clustering
  cfg = ft_checkconfig(cfg, 'unused', {'clusterstatistic', 'clusteralpha', 'clustercritval', 'clusterthreshold', 'clustertail'});
end

if any(strcmp(cfg.correctm, {'cluster' 'tfce'}))
  % these options might require a spatial neighbourhood matrix
  
  % deal with the neighbourhood of the channels/triangulation/voxels
  if isempty(cfg.connectivity)
    if isfield(cfg, 'dim') && ~isfield(cfg, 'channel') && ~isfield(cfg, 'tri')
      % input data can be reshaped into a 3D volume, use bwlabeln/spm_bwlabel rather than clusterstat
      ft_info('using connectivity of voxels in 3-D volume\n');
      cfg.connectivity = nan;
    elseif isfield(cfg, 'tri')
      % input data describes a surface along which neighbours can be defined
      ft_info('using connectivity of vertices along triangulated surface\n');
      cfg.connectivity = triangle2connectivity(cfg.tri);
      if isfield(cfg, 'insideorig')
        cfg.connectivity = cfg.connectivity(cfg.insideorig, cfg.insideorig);
      end
    elseif isfield(cfg, 'avgoverchan') && istrue(cfg.avgoverchan)
      % channel dimension has been averaged across, no sense in clustering across space
      cfg.connectivity = true(1);
    elseif isfield(cfg, 'channel')
      cfg.neighbours   = ft_getopt(cfg, 'neighbours', []);
      cfg.connectivity = channelconnectivity(cfg);
    else
      % there is no connectivity in the spatial dimension
      cfg.connectivity = false(size(dat,1));
    end
  else
    % use the specified connectivity: this is not fully robust because
    % there is no guarantee that the order of the spatial elements in the
    % data is the same as the order of the spatial elements in the
    % adjacency matrix
  end
end

% for backward compatibility and other warnings relating correcttail
if isfield(cfg,'correctp') && strcmp(cfg.correctp,'yes')
  ft_warning('cfg.correctp has been renamed to cfg.correcttail and the options have been changed')
  disp('setting cfg.correcttail to ''prob''')
  cfg.correcttail = 'prob';
  cfg = rmfield(cfg,'correctp');
elseif isfield(cfg,'correctp') && strcmp(cfg.correctp,'no')
  cfg = ft_checkconfig(cfg, 'renamed', {'correctp', 'correcttail'});
end
if strcmp(cfg.correcttail,'no') && cfg.tail==0 && cfg.alpha==0.05
  ft_warning('Doing a two-sided test without correcting p-values or alpha-level, p-values and alpha-level will reflect one-sided tests per tail. See http://bit.ly/2YQ1Hm8')
end

% for backward compatibility
if size(design,2)~=size(dat,2)
  design = transpose(design);
end

if ischar(cfg.ivar) && strcmp(cfg.ivar, 'all')
  cfg.ivar = 1:size(design,1);
end

% fetch function handle to the low-level statistics function
statfun = ft_getuserfun(cfg.statistic, 'statfun');
if isempty(statfun)
  ft_error('could not locate the appropriate statistics function');
else
  ft_info('using "%s" for the single-sample statistics\n', func2str(statfun));
end

% construct the resampled design matrix or data-shuffling matrix
ft_info('constructing randomized design\n');
resample = resampledesign(cfg, design);
Nrand = size(resample,1);

% most of the statfuns result in this warning, which is not interesting
ws = ft_warning('off', 'MATLAB:warn_r14_stucture_assignment');

if strcmp(cfg.correctm, 'cluster')
  % determine the critical value for cluster thresholding
  if strcmp(cfg.clusterthreshold, 'nonparametric_individual') || strcmp(cfg.clusterthreshold, 'nonparametric_common')
    ft_info('using a nonparametric threshold for clustering\n');
    cfg.clustercritval = [];  % this will be determined later
  elseif strcmp(cfg.clusterthreshold, 'parametric') && isempty(cfg.clustercritval)
    ft_info('computing a parametric threshold for clustering\n');
    tmpcfg = cfg; % the next line does not pass on non-standard options that a statfun might use
    % tmpcfg = keepfields(cfg, {'dim' 'dimord' 'clusteralpha' 'clustertail' 'ivar' 'uvar' 'cvar' 'wvar' 'contrastcoefs'});
    tmpcfg.computecritval = 'yes';  % explicitly request the computation of the crtitical value
    tmpcfg.computestat    = 'no';   % skip the computation of the statistic
    tmpcfg.alpha          = cfg.clusteralpha; % the statfun uses cfg.alpha most likely 
    try
      cfg.clustercritval    = getfield(statfun(tmpcfg, dat, design), 'critval');
    catch
      disp(lasterr);
      ft_error('could not determine the parametric critical value for clustering');
    end
  elseif strcmp(cfg.clusterthreshold, 'parametric') && ~isempty(cfg.clustercritval)
    ft_info('using the specified parametric threshold for clustering\n');
    cfg.clusteralpha = [];
  end
end

% compute the statistic for the observed data
ft_progress('init', cfg.feedback, 'computing statistic');
% get an estimate of the time required per evaluation of the statfun
time_pre = cputime;

try
  % the nargout function in MATLAB 6.5 and older does not work on function handles
  num = nargout(statfun);
catch
  num = 1;
end

if num==1
  % only the statistic is returned
  [statobs] = statfun(cfg, dat, design);
elseif num==2
  % both the statistic and the (updated) configuration are returned
  [statobs, cfg] = statfun(cfg, dat, design);
elseif num==3
  % both the statistic and the (updated) configuration and the (updated) data are returned
  tmpcfg = cfg;
  if strcmp(cfg. precondition, 'before'), tmpcfg.preconditionflag = 1; end
  [statobs, tmpcfg, dat]  = statfun(tmpcfg, dat, design);
  tmpcfg.preconditionflag = 0;
  cfg = tmpcfg;
end

if isstruct(statobs)
  % remember all details for later reference, continue to work with the statistic
  statfull = statobs;
  statobs  = statobs.stat;
end

% remember the statistic for later reference, continue to work with the statistic
statfull.stat = statobs;

time_eval = cputime - time_pre;
ft_info('estimated time per randomization is %.2f seconds\n', time_eval);

% pre-allocate some memory
if strcmp(cfg.correctm, 'cluster')
  statrand = zeros(size(statobs,1), size(resample,1), class(dat)); % this reduces the memory footprint, requires the user to use ft_struct2single on the input data
else
  prb_pos   = zeros(size(statobs));
  prb_neg   = zeros(size(statobs));
end

if strcmp(cfg.precondition, 'after')
  tmpcfg = cfg;
  tmpcfg.preconditionflag = 1;
  [tmpstat, tmpcfg, dat] = statfun(tmpcfg, dat, design);
end

if any(strcmp(cfg.correctm, {'tfce' 'max'}))
  % pre-allocate the memory to hold the distribution of most extreme positive (right) and negative (left) statistical values
  posdistribution = nan(1,Nrand);
  negdistribution = nan(1,Nrand);
end

% compute the statistic for the randomized data and count the outliers
for i = 1:Nrand
  ft_progress(i/Nrand, 'computing statistic %d from %d\n', i, Nrand);
  if strcmp(cfg.resampling, 'permutation')
    tmpdesign = design(:,resample(i,:));     % the columns in the design matrix are reshufled by means of permutation
    tmpdat    = dat;                        % the data itself is not shuffled
    if size(tmpdesign,1)==size(tmpdat,2)
      tmpdesign = transpose(tmpdesign);
    end
  elseif strcmp(cfg.resampling, 'bootstrap')
    tmpdesign = design;                     % the design matrix is not shuffled
    tmpdat    = dat(:,resample(i,:));        % the columns of the data are resampled by means of bootstrapping
  end
  if any(strcmp(cfg.correctm, {'cluster' 'tfce'}))
    % keep each randomization in memory for cluster postprocessing
    dum = statfun(cfg, tmpdat, tmpdesign);
    if isstruct(dum)
      statrand(:,i) = dum.stat;
    else
      statrand(:,i) = dum;
    end
  else
    % do not keep each randomization in memory, but process them on the fly
    statrand = statfun(cfg, tmpdat, tmpdesign);
    if isstruct(statrand)
      statrand = statrand.stat;
    end
   
    % the following line is for debugging
    % stat.statkeep(:,i) = statrand;
    if strcmp(cfg.correctm, 'max')
      % compare each data element with the maximum statistic
      prb_pos = prb_pos + (statobs<max(statrand(:)));
      prb_neg = prb_neg + (statobs>min(statrand(:)));
      posdistribution(i) = max(statrand(:));
      negdistribution(i) = min(statrand(:));
    else
      % compare each data element with its own statistic
      prb_pos = prb_pos + (statobs<statrand);
      prb_neg = prb_neg + (statobs>statrand);
    end
  end
end
ft_progress('close');

if strcmp(cfg.correctm, 'cluster')
  % do the cluster postprocessing
  [stat, cfg] = clusterstat(cfg, statrand, statobs);
elseif strcmp(cfg.correctm, 'tfce')
  [stat, cfg] = tfcestat(cfg, statrand, statobs);
else
  if ~isequal(cfg.numrandomization, 'all')
    % in case of random permutations (i.e., montecarlo sample, and NOT full
    % permutation), the minimum p-value should not be 0, but 1/N
    prb_pos = prb_pos + 1;
    prb_neg = prb_neg + 1;
    Nrand = Nrand + 1;
  end
  switch cfg.tail
    case 1
      clear prb_neg  % not needed any more, free some memory
      stat.prob = prb_pos./Nrand;
    case -1
      clear prb_pos  % not needed any more, free some memory
      stat.prob = prb_neg./Nrand;
    case 0
      % for each observation select the tail that corresponds with the lowest probability
      prb_neg = prb_neg./Nrand;
      prb_pos = prb_pos./Nrand;
      stat.prob = min(prb_neg, prb_pos); % this is the probability for the most unlikely tail
  end
end

% In case of a two tailed test, the type-I error rate (alpha) refers to
% both tails of the distribution, whereas the stat.prob value computed sofar
% corresponds with one tail, i.e. the probability, under the assumption of
% no effect or no difference (the null hypothesis), of obtaining a result
% equal to or more extreme than what was actually observed. The decision
% rule whether the null-hopothesis should be rejected given the observed
% probability therefore should consider alpha divided by two, to correspond
% with the probability in one of the tails (the most extreme tail). This
% is conceptually equivalent to performing a Bonferroni correction for the
% two tails.
%
% An alternative solution to distribute the alpha level over both tails is
% achieved by multiplying the probability with a factor of two, prior to
% thresholding it wich cfg.alpha.  The advantage of this solution is that
% it results in a p-value that corresponds with a parametric probability.
% Below both options are realized
if strcmp(cfg.correcttail, 'prob') && cfg.tail==0
  stat.prob = stat.prob .* 2;
  stat.prob(stat.prob>1) = 1; % clip at p=1
  % also correct the probabilities in the pos/negcluster fields
  if isfield(stat, 'posclusters')
    for i=1:length(stat.posclusters)
      stat.posclusters(i).prob = stat.posclusters(i).prob*2;
      if stat.posclusters(i).prob>1; stat.posclusters(i).prob = 1; end
    end
  end
  if isfield(stat, 'negclusters')
    for i=1:length(stat.negclusters)
      stat.negclusters(i).prob = stat.negclusters(i).prob*2;
      if stat.negclusters(i).prob>1; stat.negclusters(i).prob = 1; end
    end
  end
elseif strcmp(cfg.correcttail, 'alpha') && cfg.tail==0
  cfg.alpha = cfg.alpha / 2;
end

% compute range of confidence interval p ? 1.96(sqrt(var(p))), with var(p) = var(x/n) = p*(1-p)/N
stddev = sqrt(stat.prob.*(1-stat.prob)/Nrand);
stat.cirange = 1.96*stddev;

if isfield(stat, 'posclusters')
  for i=1:length(stat.posclusters)
    stat.posclusters(i).stddev  = sqrt(stat.posclusters(i).prob.*(1-stat.posclusters(i).prob)/Nrand);
    stat.posclusters(i).cirange =  1.96*stat.posclusters(i).stddev;
    if i==1 && stat.posclusters(i).prob<cfg.alpha && stat.posclusters(i).prob+stat.posclusters(i).cirange>=cfg.alpha
      ft_warning('FieldTrip:posCluster_exceeds_alpha', sprintf('The p-value confidence interval of positive cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
    end
  end
end
if isfield(stat, 'negclusters')
  for i=1:length(stat.negclusters)
    stat.negclusters(i).stddev  = sqrt(stat.negclusters(i).prob.*(1-stat.negclusters(i).prob)/Nrand);
    stat.negclusters(i).cirange =  1.96*stat.negclusters(i).stddev;
    if i==1 && stat.negclusters(i).prob<cfg.alpha && stat.negclusters(i).prob+stat.negclusters(i).cirange>=cfg.alpha
      ft_warning('FieldTrip:negCluster_exceeds_alpha', sprintf('The p-value confidence interval of negative cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
    end
  end
end

if ~isfield(stat, 'prob')
  ft_warning('probability was not computed');
else
  switch lower(cfg.correctm)
    case 'max'
      % the correction is implicit in the method
      ft_notice('using a maximum-statistic based method for multiple comparison correction\n');
      ft_notice('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
      stat.mask = stat.prob<=cfg.alpha;
      stat.posdistribution = posdistribution;
      stat.negdistribution = negdistribution;
    case 'tfce'
      ft_notice('using a threshold free cluster enhancement based method for multiple comparison correction\n');
      ft_notice('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
      stat.mask = stat.prob<=cfg.alpha;
    case 'cluster'
      % the correction is implicit in the method
      ft_notice('using a cluster-based method for multiple comparison correction\n');
      ft_notice('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
      stat.mask = stat.prob<=cfg.alpha;
    case 'bonferroni'
      ft_notice('performing Bonferroni correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
    case 'holm'
      % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
      ft_notice('performing Holm-Bonferroni correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
      k = find(pvals > (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'first'); % compare each significance probability against its individual threshold
      mask = (1:length(pvals))'<k;
      stat.mask = zeros(size(stat.prob));
      stat.mask(indx) = mask;
    case 'hochberg'
      % test the most significant significance probability against alpha/N, the second largest against alpha/(N-1), etc.
      ft_notice('performing Hochberg''s correction for multiple comparisons (this is *not* the Benjamini-Hochberg FDR procedure!)\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      [pvals, indx] = sort(stat.prob(:));                     % this sorts the significance probabilities from smallest to largest
      k = find(pvals <= (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'last'); % compare each significance probability against its individual threshold
      mask = (1:length(pvals))'<=k;
      stat.mask = zeros(size(stat.prob));
      stat.mask(indx) = mask;
    case 'fdr'
      ft_notice('performing FDR correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      stat.mask = fdr(stat.prob, cfg.alpha);
    otherwise
      ft_notice('not performing a correction for multiple comparisons\n');
      stat.mask = stat.prob<=cfg.alpha;
  end
end

% return the observed statistic
if ~isfield(stat, 'stat')
  stat.stat = statobs;
end

if exist('statrand', 'var')
  stat.ref = mean(statrand,2);
end

% return optional other details that were returned by the statfun
stat = copyfields(statfull, stat, fieldnames(statfull));

ft_warning(ws); % revert to original state

