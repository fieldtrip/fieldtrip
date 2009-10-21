function [stat] = sourcestatistics_randcluster(cfg,source)

% SOURCESTATISTICS_RANDCLUSTER performs statistics on the output of SOURCEANALYSIS,
% provided the randomization-method has been applied. It will be called from within the
% SOURCESTATISTICS function, when cfg.method is set to 'randcluster'. See also SOURCESTATISTICS
% for additional information.
%
% Use as
%  [stat] = sourcestatistics(cfg,source)
% where cfg is a structure containing the configuration details, and source is a structure containing
% the fields trialA and trialB, obtained with the 'randomization'-method in sourceanalysis,
%  or
%  [stat] = sourcestatistics(cfg,stat)
% where cfg is a structure containing the configuration details, and stat is the output of 
% sourcestatistics, by using the 'randcluster'-method, and using the 'intermediate'-option.
%
% For each randomization repetition, the difference volume between condition A and B is computed.
% In each volume p-values are computed for each voxel, with respect to the randomization distribution
% for that particular voxel (uncorrected for multiple comparisons). Subsequently, a clustering 
% algorithm (with a voxel connectivity of 6) locates the biggest cluster, within each randomization
% volume, with p-values smaller than a predefined value. This yields a reference distribution of
% cluster sizes, against which the size of the clusters in the observed data (by using the same cluster
% threshold) are tested. In this way, correction for multiple comparisons is achieved.
% 
% The following parameters are mandatory:
%    cfg.method           = 'randcluster'
%    cfg.parameter        = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
% The following parameters are optional:
%    cfg.comparestat      = 'difference' (default), or 'relchange' (relative change)
%    cfg.clusterthreshold = the a priori threshold at which the individual volumes will be 
%                           thresholded (default = 0.05);
%    cfg.threshold        = the p-value with respect to the reference distribution at which 
%                           an observed cluster will be considered significant. (default = 0.05)
%    cfg.intermediate     = if set to 'yes' (default), it keeps the uncorrected p-values in
%                           the output structure, so that the function's output can directly
%                           be used again as an input, which will save quite some computation time
%                           if you want to evaluate the effect of a different clusterthreshold.
%    cfg.tail             = 0, -1, or 1.
%                            if cfg.tail = 0, the alternative for a rejected null-hypothesis will be:
%                                             condition A and B are different. (default)
%                            if cfg.tail = -1,the alternative for a rejected null-hypothesis will be:
%                                             condition A < B.
%                            if cfg.tail = 1, the alternative for a rejected null-hypothesis will be:
%                                             condition A > B.                            
%    cfg.ztransform       = 'no' (default), or 'yes'. This z-transforms the tested parameter based on the
%                           randomization variance estimate for each voxel.
% The function outputs the structure stat, containing the following fields:
%     stat.cluster:      information about the clusters in the observed data
%                          nVox:         the number of voxels.
%                          prob:         the p-value against the reference distribution.
%                          shape:        the shape of the cluster, voxels belonging to the cluster
%                                         are set to 1, the rest will be 0.
%                          significance: the significance of the cluster against cfg.threshold.
%     stat.dist:         the reference distribution of maximum cluster-sizes against which the observed
%                          data has been tested.
%     stat.mask:         a volume-reconstruction with the same size as source.dim, with the significant
%                          clusters put to 1, and the rest will be 0. Can be used as a mask for displaying
%                          purposes.
%     stat.intermediate: (if specified in the configuration)
%                          prob_obs is a vector containing for each voxel the p-values (uncorrected) against
%                           the randomization-distribution for that particular voxel.
%                          prob_rand is a matrix of number of randomizations times number of voxels, containing
%                           for each randomization and voxel the p-values (uncorrected) against the randomization
%                           distribution for that particular voxel.

% FIXME this function should use parameterselection and getsubfield

% Copyright (C) 2004, Jan-Mathijs Schoffelen
%
% $Log: sourcestatistics_randcluster.m,v $
% Revision 1.15  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.14  2007/05/30 07:04:05  roboos
% use the checkdata function to avlidate the input and to convert the inside vector to indices
%
% Revision 1.13  2006/07/05 10:21:56  roboos
% updaed documentation for consistency
%
% Revision 1.12  2006/03/30 12:24:33  roboos
% Implemented private/fixinside, which facilitates consistent
% handling of source/volume data. Improved documentation. Fixed some
% bugs related to inconsistent handling of ROIs (i.e. inside/outside)
%
% Revision 1.11  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.10  2005/06/03 08:58:07  roboos
% transfer homogenous transformation matrix from input to output (if present)
% added an extra method of counting the number of dipoles/voxels in the input source structure
%
% Revision 1.9  2005/05/17 17:50:39  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.8  2005/05/12 12:17:23  roboos
% added a reshape around the source parameter to ensure that it is a row vector
% if source.pos is not available, determine the number of dipoles from source.dim
%
% Revision 1.7  2005/05/09 14:14:46  roboos
% made the assignment of the output stat fields that describe the grid (like pos, xgrid, dim ...) optional, so that it won't fail when one of these fields is not present
%
% Revision 1.6  2005/03/18 14:36:15  roboos
% fixed bug in two sided testing by replacing the relevant code by the RANDSTATPROB private subfunction
%
% Revision 1.5  2004/11/09 15:23:10  jansch
% fixed nasty bug which led to systematic underestimation of the max-cluster-
% size in the randomized volumes.
%
% Revision 1.4  2004/10/20 10:37:55  roboos
% fixed bug that occurred if avg.pow was not a row-vector
%
% Revision 1.3  2004/09/22 14:09:52  roboos
% renamed method=randomized into method=randomization, with backward compatibility
%
% Revision 1.2  2004/09/08 13:31:56  jansch
% changed some code. added version information
%

fieldtripdefs

% check if the input data is valid for this function
source = checkdata(source, 'datatype', {'source', 'volume'}, 'feedback', 'no', 'inside', 'index');

% set the defaults
if ~isfield(cfg,'threshold'),        cfg.threshold        = 0.05;         end
if ~isfield(cfg,'clusterthreshold'), cfg.clusterthreshold = 0.05;         end
if ~isfield(cfg,'tail'),             cfg.tail             = 0;            end 
if ~isfield(cfg,'intermediate'),     cfg.intermediate     = 'yes';        end
if ~isfield(cfg,'comparestat'),      cfg.comparestat      = 'difference'; end
if ~isfield(cfg, 'ztransform'),      cfg.ztransform       = 'no';         end

% this is required for backward compatibility with the old sourceanalysis
if isfield(source, 'method') && strcmp(source.method, 'randomized')
  source.method = 'randomization';
elseif isfield(source, 'method') && strcmp(source.method, 'permuted')
  source.method = 'permutation';
end

if ~isfield(cfg,'parameter')
  error('no parameter to do statistics on has been specified');
end

if ~strcmp(cfg.intermediate,'yes') && ~strcmp(source.method,'randomization'),
    error('the data are not suited to do randomization statistics on');
end

if ~isfield(source,'intermediate'),
  if isfield(source, 'pos')
    % count the number of dipole positions, which can be either on a regular or an irregular grid
    nDipole = size(source.pos,1);
  elseif isfield(source, 'dim')
    % individual dipole positions are not available after spatial normalisation to a template anatomical MRI
    nDipole = prod(source.dim);
  else
    % count the number of functional values, which should be the same as the number of dipoles
    dum = getfield(source.avgA,cfg.parameter);
    nDipole = prod(size(dum));
  end

  % get the data to work on in nice arrays
  nTrial  = size(source.trialA,2);
  strialA = zeros(nTrial,nDipole);
  strialB = zeros(nTrial,nDipole);
  for j = 1:nTrial
   % reshape to ensure that the functional parameter for all voxels is arranged in a row-vector
   strialA(j,:) = reshape(getfield(source.trialA,{j},cfg.parameter), [1 nDipole]);
   strialB(j,:) = reshape(getfield(source.trialB,{j},cfg.parameter), [1 nDipole]);
  end
  % reshape to ensure that the functional parameter for all voxels is arranged in a row-vector
  savgA  = reshape(getfield(source.avgA,cfg.parameter), [1 nDipole]);
  savgB  = reshape(getfield(source.avgB,cfg.parameter), [1 nDipole]);
  % transform to z-scores if requested
  if strcmp(cfg.ztransform,'yes')
    inside = source.inside;
    meanA = repmat(mean(strialA(:,inside)),size(strialA,1),1);
    meanB = repmat(mean(strialB(:,inside)),size(strialB,1),1);
    stdA = repmat(std(strialA(:,inside)),size(strialA,1),1);
    stdB = repmat(std(strialB(:,inside)),size(strialB,1),1);
    strialA(:,inside) = (strialA(:,inside) - meanA) ./ stdA; 
    strialB(:,inside) = (strialB(:,inside) - meanB) ./ stdB;
    savgA(1,inside) = (savgA(1,inside) - meanA(1,:)) ./ stdA(1,:);
    savgB(1,inside) = (savgB(1,inside) - meanB(1,:)) ./ stdB(1,:);
  end
  % normalize by condition B if requested
  if strcmp(cfg.comparestat,'difference'),
     randobs     = strialA-strialB;
     realobs     = savgA-savgB;
  elseif strcmp(cfg.comparestat,'relchange'),
     randobs     = (strialA-strialB)./strialB;
     realobs     = (savgA-savgB)./savgB;
  end;
else
  % only the number of trials needs to be determined, the preprocessing up to
  % the per-voxel-probability computation already have been performed before
  nTrial = size(source.trialA,2);
end

inside  = source.inside;
outside = source.outside;
c     = zeros(nTrial,1);
maxC  = zeros(nTrial,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the maximum-cluster size on the thresholded p-values
%for each randomization
%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:nTrial
  fprintf('trial %d of %d\n',k,nTrial);

  if ~isfield(source,'intermediate'),
     % compute p-values
     intermediate.prob_rand(k,inside)  = randstatprob(randobs(:,inside)', randobs(k,inside)', cfg.tail, 0)';
     intermediate.prob_rand(k,outside) = nan;
  else
     % reuse the previously computed p-values
     intermediate.prob_rand = source.intermediate.prob_rand;
  end

  dum = intermediate.prob_rand(k,:);     % take the probability over all voxels
  dum = find(dum<cfg.clusterthreshold);  % threshold it and find the clusters
  onoff = zeros(source.dim);
  onoff(dum) = 1;
  if ~isempty(dum),
      remember = findconnect(onoff);
      %remember = findcluster(onoff);
      numC    = length(remember);
      maxc    = 0;
      for j = 1:numC
        maxc = max([maxc length(find(remember{j}))]);
      end
      maxC(k) = maxc;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the clusters of thresholded p-values
%in the observed data
%%%%%%%%%%%%%%%%%%%%%

if ~isfield(source,'intermediate')
  intermediate.prob_obs(:,inside)  = randstatprob(randobs(:,inside)', realobs(:,inside)', cfg.tail, 0)';
  intermediate.prob_obs(:,outside) = nan;
else
  intermediate.prob_obs = source.intermediate.prob_obs;
end

dum = intermediate.prob_obs;           % take the probability over all voxels
dum = find(dum<cfg.clusterthreshold);  % threshold it and find the clusters
onoff = zeros(source.dim);
onoff(dum) = 1;
if ~isempty(dum),
    remember = findconnect(onoff);
    %remember = findcluster(onoff);
    numC     = length(remember);
    mask     = logical(zeros(size(onoff)));
    for k = 1:numC
     cluster(k).nVox = length(find(remember{k})); 
     cluster(k).prob = (length(find(maxC>cluster(k).nVox))+1)./(nTrial);
     cluster(k).shape= remember{k};
     cluster(k).significance  = logical(cluster(k).prob < cfg.threshold);
     if cluster(k).significance
        mask = (mask | cluster(k).shape);
     end 
    end
    [srt,indx] = sort([cluster(:).nVox]);
    cluster=cluster(indx);
end

%%%%%%%%%%%%%%%%%%%%
%collect the results
%%%%%%%%%%%%%%%%%%%%

try, stat.dim         = source.dim;       end
try, stat.xgrid       = source.xgrid;     end
try, stat.ygrid       = source.ygrid;     end
try, stat.zgrid       = source.zgrid;     end
try, stat.inside      = source.inside;    end
try, stat.outside     = source.outside;   end
try, stat.pos         = source.pos;       end
try, stat.transform   = source.transform; end

stat.cluster          = cluster;
stat.dist             = maxC;
stat.mask             = mask;
if strcmp(cfg.intermediate,'yes'), stat.intermediate = intermediate; end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: sourcestatistics_randcluster.m,v 1.15 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = source.cfg ; end
% remember the exact configuration details in the output 
stat.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction for finding connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function remember = findconnect(onoff)

seg = bwlabeln(onoff, 6);
for i=1:max(seg(:));
  remember{i} = (seg==i);
end
