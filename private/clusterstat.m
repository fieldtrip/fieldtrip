function [stat, cfg] = clusterstat(cfg, statrnd, statobs, varargin)

% SUBFUNCTION for computing cluster statistic for N-D volumetric source data
% or for channel-freq-time data
%
% This function uses
%   cfg.dim
%   cfg.inside (only for source data)
%   cfg.tail = -1, 0, 1
%   cfg.multivariate = no, yes
%   cfg.orderedstats = no, yes
%   cfg.clusterstatistic = max, maxsize, maxsum, wcm
%   cfg.clusterthreshold = parametric, nonparametric_individual, nonparametric_common
%   cfg.clusteralpha
%   cfg.clustercritval
%   cfg.wcm_weight
%   cfg.feedback

% Copyright (C) 2005-2007, Robert Oostenveld
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

% set the defaults
cfg.orderedstats = ft_getopt(cfg, 'orderedstats', 'no');
cfg.multivariate = ft_getopt(cfg, 'multivariate', 'no');
cfg.minnbchan    = ft_getopt(cfg, 'minnbchan',    0);

if cfg.tail~=cfg.clustertail
  error('cfg.tail and cfg.clustertail should be identical')
end

% get conncevitiy matrix for the spatially neighbouring elements
channeighbstructmat = full(ft_getopt(cfg, 'connectivity', false));

needpos = cfg.tail==0 || cfg.tail== 1;
needneg = cfg.tail==0 || cfg.tail==-1;
Nsample = size(statrnd,1);
Nrand   = size(statrnd,2);

prb_pos    = ones(Nsample,     1);
prb_neg    = ones(Nsample,     1);
postailrnd = false(Nsample,Nrand);  % this holds the thresholded values
negtailrnd = false(Nsample,Nrand);  % this holds the thresholded values
Nobspos = 0;                        % number of positive clusters in observed data
Nobsneg = 0;                        % number of negative clusters in observed data

if strcmp(cfg.clusterthreshold, 'parametric')
  % threshold based on the critical value from parametric distribution
  siz = size(cfg.clustercritval);
  if all(siz==1) && cfg.clustertail==0
    %  it only specifies one critical value, assume that the left and right tail are symmetric around zero
    negtailcritval = -cfg.clustercritval;
    postailcritval =  cfg.clustercritval;
  elseif all(siz==1) && cfg.clustertail==-1
    % it only specifies one critical value corresponding to the left tail
    negtailcritval = cfg.clustercritval;
    postailcritval = +inf * ones(size(negtailcritval));
  elseif all(siz==1) && cfg.clustertail==1
    % it only specifies one critical value corresponding to the right tail
    postailcritval =  cfg.clustercritval;
    negtailcritval = -inf * ones(size(postailcritval));
  elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==0
    %  it specifies a single critical value for each sample, assume that the left and right tail are symmetric around zero
    negtailcritval = -cfg.clustercritval;
    postailcritval =  cfg.clustercritval;
  elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==-1
    % it specifies a critical value for the left tail
    % which is different for each sample (samples have a different df)
    negtailcritval = cfg.clustercritval;
    postailcritval = +inf * ones(size(negtailcritval));
  elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==1
    % it specifies a critical value for the right tail
    % which is different for each sample (samples have a different df)
    postailcritval = cfg.clustercritval;
    negtailcritval = +inf * ones(size(postailcritval));
  elseif siz(1)==Nsample && siz(2)==2 && cfg.clustertail==0
    % it specifies a critical value for the left and for the right tail of the distribution
    % which is different for each sample (samples have a different df)
    negtailcritval = cfg.clustercritval(:,1);
    postailcritval = cfg.clustercritval(:,2);
  elseif prod(siz)==2 && cfg.clustertail==0
    % it specifies a critical value for the left and for the right tail of the distribution
    % which is the same for each sample (samples have the same df)
    negtailcritval = cfg.clustercritval(1);
    postailcritval = cfg.clustercritval(2);
  else
    error('cannot make sense out of the specified parametric critical values');
  end
  
elseif strcmp(cfg.clusterthreshold, 'nonparametric_individual')
  % threshold based on bootstrap using all other randomizations
  % each voxel will get an individual threshold
  [srt, ind] = sort(statrnd,2);
  if cfg.clustertail==0
    % both tails are needed
    negtailcritval = srt(:,round((  cfg.clusteralpha/2)*size(statrnd,2)));
    postailcritval = srt(:,round((1-cfg.clusteralpha/2)*size(statrnd,2)));
  elseif cfg.clustertail==1
    % only positive tail is needed
    postailcritval = srt(:,round((1-cfg.clusteralpha)*size(statrnd,2)));
    negtailcritval = -inf * ones(size(postailcritval));
  elseif cfg.clustertail==-1
    % only negative tail is needed
    negtailcritval = srt(:,round((  cfg.clusteralpha)*size(statrnd,2)));
    postailcritval = +inf * ones(size(negtailcritval));
  end
  
elseif strcmp(cfg.clusterthreshold, 'nonparametric_common')
  % threshold based on bootstrap using all other randomizations
  % all voxels will get a common threshold
  [srt, ind] = sort(statrnd(:));
  if cfg.clustertail==0
    % both tails are needed
    negtailcritval = srt(round((  cfg.clusteralpha/2)*prod(size(statrnd))));
    postailcritval = srt(round((1-cfg.clusteralpha/2)*prod(size(statrnd))));
  elseif cfg.clustertail==1
    % only positive tail is needed
    postailcritval = srt(round((1-cfg.clusteralpha)*prod(size(statrnd))));
    negtailcritval = -inf * ones(size(postailcritval));
  elseif cfg.clustertail==-1
    % only negative tail is needed
    negtailcritval = srt(round((  cfg.clusteralpha)*prod(size(statrnd))));
    postailcritval = +inf * ones(size(negtailcritval));
  end
  
else
  error('no valid threshold for clustering was given')
end % determine clusterthreshold

% these should be scalars or column vectors
negtailcritval = negtailcritval(:);
postailcritval = postailcritval(:);

% remember the critical values
cfg.clustercritval = [negtailcritval postailcritval];

% test whether the observed and the random statistics exceed the threshold
postailobs = (statobs >= postailcritval);
negtailobs = (statobs <= negtailcritval);
for i=1:Nrand
  postailrnd(:,i) = (statrnd(:,i) >= postailcritval);
  negtailrnd(:,i) = (statrnd(:,i) <= negtailcritval);
end

% first do the clustering on the observed data
if needpos,
  
  if ~isfinite(channeighbstructmat)
    % this pertains to data for which the spatial dimension can be reshaped
    % into 3D, i.e. when it is described on an ordered set of positions on a 3D-grid
    if isfield(cfg, 'origdim'), 
      cfg.dim = cfg.origdim; 
    end %this snippet is to support correct clustering of N-dimensional data, not fully tested yet
    tmp = zeros(cfg.dim);
    tmp(cfg.inside) = postailobs;
    
    numdims = length(cfg.dim);
    if numdims == 2 || numdims == 3 % if 2D or 3D data
      ft_hastoolbox('spm8',1);
      [posclusobs, posnum] = spm_bwlabel(tmp, 2*numdims); % use spm_bwlabel for 2D/3D data to avoid usage of image toolbox
    else
      posclusobs = bwlabeln(tmp, conndef(length(cfg.dim),'min')); % spm_bwlabel yet (feb 2011) supports only 2D/3D data
    end
    posclusobs = posclusobs(cfg.inside);
  else
    if 0
      posclusobs = findcluster(reshape(postailobs, [cfg.dim,1]),cfg.chancmbneighbstructmat,cfg.chancmbneighbselmat,cfg.minnbchan);
    else
      posclusobs = findcluster(reshape(postailobs, [cfg.dim,1]),channeighbstructmat,cfg.minnbchan);
    end
    posclusobs = posclusobs(:);
  end
  Nobspos = max(posclusobs(:)); % number of clusters exceeding the threshold
  fprintf('found %d positive clusters in observed data\n', Nobspos);

end
if needneg,
  
  if ~isfinite(channeighbstructmat)
    % this pertains to data for which the spatial dimension can be reshaped
    % into 3D, i.e. when it is described on an ordered set of positions on a 3D-grid
    
    tmp = zeros(cfg.dim);
    tmp(cfg.inside) = negtailobs;
      
    numdims = length(cfg.dim);
    if numdims == 2 || numdims == 3 % if 2D or 3D data
      ft_hastoolbox('spm8',1);
      [negclusobs, negnum] = spm_bwlabel(tmp, 2*numdims); % use spm_bwlabel for 2D/3D data to avoid usage of image toolbox
    else
      negclusobs = bwlabeln(tmp, conndef(length(cfg.dim),'min')); % spm_bwlabel yet (feb 2011) supports only 2D/3D data
    end
    negclusobs = negclusobs(cfg.inside);
  else
    if 0
      negclusobs = findcluster(reshape(negtailobs, [cfg.dim,1]),cfg.chancmbneighbstructmat,cfg.chancmbneighbselmat,cfg.minnbchan);
    else
      negclusobs = findcluster(reshape(negtailobs, [cfg.dim,1]),channeighbstructmat,cfg.minnbchan);
    end
    negclusobs = negclusobs(:);
  end
  Nobsneg = max(negclusobs(:));
  fprintf('found %d negative clusters in observed data\n', Nobsneg);

end

stat = [];
stat.stat = statobs;

% catch situation where no clustering of the random data is needed
if (Nobspos+Nobsneg)==0
  warning('no clusters were found in the observed data');
  stat.prob = ones(Nsample, 1);
  return
end

% allocate space to hold the randomization distributions of the cluster statistic
if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
  fprintf('allocating space for a %d-multivariate distribution of the positive clusters\n', Nobspos);
  fprintf('allocating space for a %d-multivariate distribution of the negative clusters\n', Nobsneg);
  posdistribution = zeros(Nobspos,Nrand); % this holds the multivariate randomization distribution of the positive cluster statistics
  negdistribution = zeros(Nobsneg,Nrand); % this holds the multivariate randomization distribution of the negative cluster statistics
else
  posdistribution = zeros(1,Nrand);       % this holds the statistic of the largest positive cluster in each randomization
  negdistribution = zeros(1,Nrand);       % this holds the statistic of the largest negative cluster in each randomization
end

% do the clustering on the randomized data
ft_progress('init', cfg.feedback, 'computing clusters in randomization');
for i=1:Nrand
  ft_progress(i/Nrand, 'computing clusters in randomization %d from %d\n', i, Nrand);
  if needpos,
    if ~isfinite(channeighbstructmat)
      tmp = zeros(cfg.dim);
      tmp(cfg.inside) = postailrnd(:,i);
      
      numdims = length(cfg.dim);
      if numdims == 2 || numdims == 3 % if 2D or 3D data
        [posclusrnd, posrndnum] = spm_bwlabel(tmp, 2*numdims); % use spm_bwlabel for 2D/3D data to avoid usage of image toolbox
      else
        posclusrnd = bwlabeln(tmp, conndef(length(cfg.dim),'min')); % spm_bwlabel yet (feb 2011) supports only 2D/3D data
      end
      posclusrnd = posclusrnd(cfg.inside);
    else
      if 0
        posclusrnd = findcluster(reshape(postailrnd(:,i), [cfg.dim,1]),cfg.chancmbneighbstructmat,cfg.chancmbneighbselmat,cfg.minnbchan);
      else
        posclusrnd = findcluster(reshape(postailrnd(:,i), [cfg.dim,1]),channeighbstructmat,cfg.minnbchan);
      end
      posclusrnd = posclusrnd(:);
    end
    Nrndpos = max(posclusrnd(:));  % number of clusters exceeding the threshold
    stat    = zeros(1,Nrndpos); % this will hold the statistic for each cluster
    % fprintf('found %d positive clusters in this randomization\n', Nrndpos);
    for j = 1:Nrndpos
      if strcmp(cfg.clusterstatistic, 'max'),
        stat(j) = max(statrnd(find(posclusrnd==j),i));
      elseif strcmp(cfg.clusterstatistic, 'maxsize'),
        stat(j) = length(find(posclusrnd==j));
      elseif strcmp(cfg.clusterstatistic, 'maxsum'),
        stat(j) = sum(statrnd(find(posclusrnd==j),i));
      elseif strcmp(cfg.clusterstatistic, 'wcm'),
        stat(j) = sum((statrnd(find(posclusrnd==j),i)-postailcritval).^cfg.wcm_weight);
      else
        error('unknown clusterstatistic');
      end
    end % for 1:Nrdnpos
    if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
      stat = sort(stat, 'descend');             % sort them from most positive to most negative
      if Nrndpos>Nobspos
        posdistribution(:,i) = stat(1:Nobspos); % remember the largest N clusters
      else
        posdistribution(1:Nrndpos,i) = stat;    % remember the largest N clusters
      end
    else
      % univariate -> remember the most extreme cluster
      if ~isempty(stat), posdistribution(i) = max(stat); end
    end
  end % needpos
  if needneg,
    if ~isfinite(channeighbstructmat)
      
      tmp = zeros(cfg.dim);
      tmp(cfg.inside) = negtailrnd(:,i);
        
      numdims = length(cfg.dim);
      if numdims == 2 || numdims == 3 % if 2D or 3D data
        ft_hastoolbox('spm8',1);
        [negclusrnd, negrndnum] = spm_bwlabel(tmp, 2*numdims); % use spm_bwlabel for 2D/3D to avoid usage of image toolbox
      else
        negclusrnd = bwlabeln(tmp, conndef(length(cfg.dim),'min')); % spm_bwlabel yet (feb 2011) supports only 2D/3D data
      end
      negclusrnd = negclusrnd(cfg.inside);
    else
      if  0
        negclusrnd = findcluster(reshape(negtailrnd(:,i), [cfg.dim,1]),cfg.chancmbneighbstructmat,cfg.chancmbneighbselmat,cfg.minnbchan);
      else
        negclusrnd = findcluster(reshape(negtailrnd(:,i), [cfg.dim,1]),channeighbstructmat,cfg.minnbchan);
      end
      negclusrnd = negclusrnd(:);
    end
    Nrndneg = max(negclusrnd(:));  % number of clusters exceeding the threshold
    stat    = zeros(1,Nrndneg); % this will hold the statistic for each cluster
    % fprintf('found %d negative clusters in this randomization\n', Nrndneg);
    for j = 1:Nrndneg
      if strcmp(cfg.clusterstatistic, 'max'),
        stat(j) = min(statrnd(find(negclusrnd==j),i));
      elseif strcmp(cfg.clusterstatistic, 'maxsize'),
        stat(j) = -length(find(negclusrnd==j)); % encode the size of a negative cluster as a negative value
      elseif strcmp(cfg.clusterstatistic, 'maxsum'),
        stat(j) = sum(statrnd(find(negclusrnd==j),i));
      elseif strcmp(cfg.clusterstatistic, 'wcm'),
        stat(j) = -sum((abs(statrnd(find(negclusrnd==j),i)-negtailcritval)).^cfg.wcm_weight); % encoded as a negative value
      else
        error('unknown clusterstatistic');
      end
    end % for 1:Nrndneg
    if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
      stat = sort(stat, 'ascend');              % sort them from most negative to most positive
      if Nrndneg>Nobsneg
        negdistribution(:,i) = stat(1:Nobsneg); % remember the most extreme clusters, i.e. the most negative
      else
        negdistribution(1:Nrndneg,i) = stat;    % remember the most extreme clusters, i.e. the most negative
      end
    else
      % univariate -> remember the most extreme cluster, which is the most negative
      if ~isempty(stat), negdistribution(i) = min(stat); end
    end
  end % needneg
end % for 1:Nrand
ft_progress('close');

% compare the values for the observed clusters with the randomization distribution
if needpos,
  posclusters = [];
  stat = zeros(1,Nobspos);
  for j = 1:Nobspos
    if strcmp(cfg.clusterstatistic, 'max'),
      stat(j) = max(statobs(find(posclusobs==j)));
    elseif strcmp(cfg.clusterstatistic, 'maxsize'),
      stat(j) = length(find(posclusobs==j));
    elseif strcmp(cfg.clusterstatistic, 'maxsum'),
      stat(j) = sum(statobs(find(posclusobs==j)));
    elseif strcmp(cfg.clusterstatistic, 'wcm'),
      stat(j) = sum((statobs(find(posclusobs==j))-postailcritval).^cfg.wcm_weight);
    else
      error('unknown clusterstatistic');
    end
  end
  % sort the clusters based on their statistical value
  [stat, indx] = sort(stat,'descend');
  % reorder the cluster indices in the data
  tmp = zeros(size(posclusobs));
  for j=1:Nobspos
    tmp(find(posclusobs==indx(j))) = j;
  end
  posclusobs = tmp;
  if strcmp(cfg.multivariate, 'yes')
    % estimate the probability of the mutivariate tail, i.e. one p-value for all clusters
    prob = 0;
    for i=1:Nrand
      % compare all clusters simultaneosuly
      prob = prob + any(posdistribution(:,i)>stat(:));
    end
    if isequal(cfg.numrandomization, 'all')
      prob = prob/Nrand;
    else % the minimum possible p-value should not be 0, but 1/N
      prob = (prob + 1)/(Nrand + 1);
    end
    for j = 1:Nobspos
      % collect a summary of the cluster properties
      posclusters(j).prob = prob;
      posclusters(j).clusterstat = stat(j);
    end
    % collect the probabilities in one large array
    prb_pos(find(posclusobs~=0)) = prob;
  elseif strcmp(cfg.orderedstats, 'yes')
    % compare the Nth ovbserved cluster against the randomization distribution of the Nth cluster
    prob = zeros(1,Nobspos);
    for j = 1:Nobspos
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(posdistribution(j,:)>stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(posdistribution(j,:)>stat(j)) + 1)/(Nrand + 1);
      end
      % collect a summary of the cluster properties
      posclusters(j).prob = prob(j);
      posclusters(j).clusterstat = stat(j);
      % collect the probabilities in one large array
      prb_pos(find(posclusobs==j)) = prob(j);
    end
  else
    % univariate -> each cluster has it's own probability
    prob = zeros(1,Nobspos);
    for j = 1:Nobspos
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(posdistribution>stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(posdistribution>stat(j)) + 1)/(Nrand + 1);
      end
      % collect a summary of the cluster properties
      posclusters(j).prob = prob(j);
      posclusters(j).clusterstat = stat(j);
      % collect the probabilities in one large array
      prb_pos(find(posclusobs==j)) = prob(j);
    end
  end
end

if needneg,
  negclusters = [];
  stat = zeros(1,Nobsneg);
  for j = 1:Nobsneg
    if strcmp(cfg.clusterstatistic, 'max'),
      stat(j) = min(statobs(find(negclusobs==j)));
    elseif strcmp(cfg.clusterstatistic, 'maxsize'),
      stat(j) = -length(find(negclusobs==j)); % encode the size of a negative cluster as a negative value
    elseif strcmp(cfg.clusterstatistic, 'maxsum'),
      stat(j) = sum(statobs(find(negclusobs==j)));
    elseif strcmp(cfg.clusterstatistic, 'wcm'),
      stat(j) = -sum((abs(statobs(find(negclusobs==j))-negtailcritval)).^cfg.wcm_weight); % encoded as a negative value
    else
      error('unknown clusterstatistic');
    end
  end
  % sort the clusters based on their statistical value
  [stat, indx] = sort(stat,'ascend');
  % reorder the cluster indices in the observed data
  tmp = zeros(size(negclusobs));
  for j=1:Nobsneg
    tmp(find(negclusobs==indx(j))) = j;
  end
  negclusobs = tmp;
  if strcmp(cfg.multivariate, 'yes')
    % estimate the probability of the mutivariate tail, i.e. one p-value for all clusters
    prob = 0;
    for i=1:Nrand
      % compare all clusters simultaneosuly
      prob = prob + any(negdistribution(:,i)<stat(:));
    end
    if isequal(cfg.numrandomization, 'all')
      prob = prob/Nrand;
    else % the minimum possible p-value should not be 0, but 1/N
      prob = (prob + 1)/(Nrand + 1);
    end
    for j = 1:Nobsneg
      % collect a summary of the cluster properties
      negclusters(j).prob = prob;
      negclusters(j).clusterstat = stat(j);
    end
    % collect the probabilities in one large array
    prb_neg(find(negclusobs~=0)) = prob;
  elseif strcmp(cfg.orderedstats, 'yes')
    % compare the Nth ovbserved cluster against the randomization distribution of the Nth cluster
    prob = zeros(1,Nobsneg);
    for j = 1:Nobsneg
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(negdistribution(j,:)<stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(negdistribution(j,:)<stat(j)) + 1)/(Nrand + 1);
      end
      % collect a summary of the cluster properties
      negclusters(j).prob = prob(j);
      negclusters(j).clusterstat = stat(j);
      % collect the probabilities in one large array
      prb_neg(find(negclusobs==j)) = prob(j);
    end
  else
    % univariate -> each cluster has it's own probability
    prob = zeros(1,Nobsneg);
    for j = 1:Nobsneg
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(negdistribution<stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(negdistribution<stat(j)) + 1)/(Nrand + 1);
      end
      % collect a summary of the cluster properties
      negclusters(j).prob = prob(j);
      negclusters(j).clusterstat = stat(j);
      % collect the probabilities in one large array
      prb_neg(find(negclusobs==j)) = prob(j);
    end
  end
end

if cfg.tail==0
  % consider both tails
  prob = min(prb_neg, prb_pos); % this is the probability for the most unlikely tail
elseif cfg.tail==1
  % only consider the positive tail
  prob = prb_pos;
elseif cfg.tail==-1
  % only consider the negative tail
  prob = prb_neg;
end

% collect the remaining details in the output structure
stat.prob = prob;
if needpos,
  stat.posclusters         = posclusters;
  stat.posclusterslabelmat = posclusobs;
  stat.posdistribution     = posdistribution;
end
if needneg,
  stat.negclusters         = negclusters;
  stat.negclusterslabelmat = negclusobs;
  stat.negdistribution     = negdistribution;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION MAKECHANNEIGHBSTRUCTMAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channeighbstructmat = makechanneighbstructmat(cfg);

% MAKECHANNEIGHBSTRUCTMAT makes the makes the matrix containing the channel
% neighbourhood structure.

% because clusterstat has no access to the actual data (containing data.label), this workaround is required
% cfg.neighbours is cleared here because it is not done where avgoverchan is effectuated (it should actually be changed there)
if strcmp(cfg.avgoverchan, 'no')
  nchan=length(cfg.channel);
elseif strcmp(cfg.avgoverchan, 'yes')
  nchan = 1;
  cfg.neighbours = [];
end
channeighbstructmat = false(nchan,nchan);
for chan=1:length(cfg.neighbours)
  [seld] = match_str(cfg.channel, cfg.neighbours(chan).label);
  [seln] = match_str(cfg.channel, cfg.neighbours(chan).neighblabel);
  if isempty(seld)
    % this channel was not present in the data
    continue;
  else
    % add the neighbours of this channel to the matrix
    channeighbstructmat(seld, seln) = true;
  end
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION MAKECHANCMBNEIGHBMATS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [chancmbneighbstructmat, chancmbneighbselmat] = makechancmbneighbmats(channeighbstructmat,labelcmb,label,orderedchancmbs,original);

% compute CHANCMBNEIGHBSTRUCTMAT and CHANCMBNEIGHBSELMAT, which will be
% used for clustering channel combinations.

nchan=length(label);
nchansqr=nchan^2;

% Construct an array labelcmbnrs from labelcmb by replacing the label pairs
% in labelcmb by numbers that correspond to the order of the labels in label.
nchancmb = size(labelcmb,1);
labelcmbnrs=zeros(nchancmb,2);
for chanindx=1:nchan
  [chansel1] = match_str(labelcmb(:,1),label(chanindx));
  labelcmbnrs(chansel1,1)=chanindx;
  [chansel2] = match_str(labelcmb(:,2),label(chanindx));
  labelcmbnrs(chansel2,2)=chanindx;
end;
% Calculate the row and column indices (which are identical) of the
% channel combinations that are present in the data.
chancmbindcs=zeros(nchancmb,1);
for indx=1:nchancmb
  chancmbindcs(indx)=(labelcmbnrs(indx,1)-1)*nchan + labelcmbnrs(indx,2);
end;

% put all elements on the diagonal of CHANNEIGHBSTRUCTMAT equal to one
channeighbstructmat=channeighbstructmat | logical(eye(nchan));

% Compute CHANCMBNEIGHBSTRUCTMAT
% First compute the complete CHANCMBNEIGHBSTRUCTMAT (containing all
% ORDERED channel combinations that can be formed with all channels present in
% the data) and later select and reorder the channel combinations actually
% present in data). In the complete CHANCMBNEIGHBSTRUCTMAT, the row and the
% column pairs are ordered lexicographically.

chancmbneighbstructmat = false(nchansqr);
if original
  % two channel pairs are neighbours if their first elements are
  % neighbours
  chancmbneighbstructmat = logical(kron(channeighbstructmat,ones(nchan)));
  
  % or if their second elements are neighbours
  chancmbneighbstructmat = chancmbneighbstructmat | logical(kron(ones(nchan),channeighbstructmat));
else    %  version that consumes less memory
  for chanindx=1:nchan
    % two channel pairs are neighbours if their first elements are neighbours
    % or if their second elements are neighbours
    chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
      logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) | ...
      logical(kron(ones(nchan,1),channeighbstructmat));
  end;
end;

if ~orderedchancmbs
  if original
    % or if the first element of the row channel pair is a neighbour of the
    % second element of the column channel pair
    chancmbneighbstructmat = chancmbneighbstructmat | logical(repmat(kron(channeighbstructmat,ones(nchan,1)), [1 nchan]));
    
    % or if the first element of the column channel pair is a neighbour of the
    % second element of the row channel pair
    chancmbneighbstructmat = chancmbneighbstructmat | logical(repmat(kron(channeighbstructmat,ones(1,nchan)),[nchan 1]));
  else
    for chanindx=1:nchan
      % two channel pairs are neighbours if
      % the first element of the row channel pair is a neighbour of the
      % second element of the column channel pair
      % or if the first element of the column channel pair is a neighbour of the
      % second element of the row channel pair
      chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
        chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) | ...
        logical(kron(channeighbstructmat,ones(nchan,1))) | ...
        logical(repmat(kron(channeighbstructmat(:,chanindx),ones(1,nchan)),[nchan 1]));
    end;
  end;
end;

% reorder and select the entries in chancmbneighbstructmat such that they correspond to labelcmb.
chancmbneighbstructmat = sparse(chancmbneighbstructmat);
chancmbneighbstructmat = chancmbneighbstructmat(chancmbindcs,chancmbindcs);

% compute CHANCMBNEIGHBSELMAT
% CHANCMBNEIGHBSELMAT identifies so-called parallel pairs. A channel pair
% is parallel if (a) all four sensors are different and (b) all elements (of the first
% and the second pair) are neighbours of an element of the other pair.
% if orderedpairs is true, then condition (b) is as follows: the first
% elements of the two pairs are neighbours, and the second elements of the
% two pairs are neighbours.

% put all elements on the diagonal of CHANNEIGHBSTRUCTMAT equal to zero
channeighbstructmat = logical(channeighbstructmat.*(ones(nchan)-diag(ones(nchan,1))));

chancmbneighbselmat = false(nchansqr);
if orderedchancmbs
  if original
    % if the first element of the row pair is a neighbour of the
    % first element of the column pair
    chancmbneighbselmat = logical(kron(channeighbstructmat,ones(nchan)));
    % and the second element of the row pair is a neighbour of the
    % second element of the column pair
    chancmbneighbselmat = chancmbneighbselmat & logical(kron(ones(nchan),channeighbstructmat));
  else    %  version that consumes less memory
    for chanindx=1:nchan
      % if the first element of the row pair is a neighbour of the
      % first element of the column pair
      % and the second element of the row pair is a neighbour of the
      % second element of the column pair
      chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
        logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) & logical(kron(ones(nchan,1),channeighbstructmat));
    end;
  end;
else  % unordered channel combinations
  if original
    % if the first element of the row pair is a neighbour of one of the
    % two elements of the column pair
    chancmbneighbselmat = logical(kron(channeighbstructmat,ones(nchan))) | logical(repmat(kron(channeighbstructmat,ones(nchan,1)), [1 nchan]));
    % and the second element of the row pair is a neighbour of one of the
    % two elements of the column pair
    chancmbneighbselmat = chancmbneighbselmat & (logical(kron(ones(nchan),channeighbstructmat)) | ...
      logical(repmat(kron(channeighbstructmat,ones(1,nchan)), [nchan 1])));
  else    %  version that consumes less memory
    for chanindx=1:nchan
      % if the first element of the row pair is a neighbour of one of the
      % two elements of the column pair
      % and the second element of the row pair is a neighbour of one of the
      % two elements of the column pair
      chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
        (logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) | logical(kron(channeighbstructmat,ones(nchan,1)))) ...
        & (logical(kron(ones(nchan,1),channeighbstructmat)) | ...
        logical(repmat(kron(channeighbstructmat(:,chanindx),ones(1,nchan)), [nchan 1])));
    end;
  end;
end;

if original
  % remove all pairs of channel combinations that have one channel in common.
  % common channel in the first elements of the two pairs.
  chancmbneighbselmat = chancmbneighbselmat & kron(~eye(nchan),ones(nchan));
  
  % common channel in the second elements of the two pairs.
  chancmbneighbselmat = chancmbneighbselmat & kron(ones(nchan),~eye(nchan));
  
  % common channel in the first element of the row pair and the second
  % element of the column pair.
  tempselmat=logical(zeros(nchansqr));
  tempselmat(:)= ~repmat([repmat([ones(nchan,1); zeros(nchansqr,1)],[(nchan-1) 1]); ones(nchan,1)],[nchan 1]);
  chancmbneighbselmat = chancmbneighbselmat & tempselmat;
  
  % common channel in the second element of the row pair and the first
  % element of the column pair.
  chancmbneighbselmat = chancmbneighbselmat & tempselmat';
else
  noteye=~eye(nchan);
  tempselmat=logical(zeros(nchansqr,nchan));
  tempselmat(:)= ~[repmat([ones(nchan,1); zeros(nchansqr,1)],[(nchan-1) 1]); ones(nchan,1)];
  for chanindx=1:nchan
    % remove all pairs of channel combinations that have one channel in common.
    % common channel in the first elements of the two pairs.
    % common channel in the second elements of the two pairs.
    % common channel in the first element of the row pair and the second
    % element of the column pair.
    chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
      chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) ...
      & logical(kron(noteye(:,chanindx),ones(nchan))) ...
      & logical(kron(ones(nchan,1),noteye)) ...
      & tempselmat;
  end;
  for chanindx=1:nchan
    % remove all pairs of channel combinations that have one
    % common channel in the second element of the row pair and the first
    % element of the column pair.
    chancmbneighbselmat(((chanindx-1)*nchan + 1):(chanindx*nchan),:) = ...
      chancmbneighbselmat(((chanindx-1)*nchan + 1):(chanindx*nchan),:) & tempselmat';
  end;
end;

% reorder and select the entries in chancmbneighbselmat such that they correspond to labelcmb.
chancmbneighbselmat = sparse(chancmbneighbselmat);
chancmbneighbselmat = chancmbneighbselmat(chancmbindcs,chancmbindcs);

% put all elements below and on the diagonal equal to zero
nchancmbindcs = length(chancmbindcs);
for chancmbindx=1:nchancmbindcs
  selvec=[true(chancmbindx-1,1);false(nchancmbindcs-chancmbindx+1,1)];
  chancmbneighbstructmat(:,chancmbindx) = chancmbneighbstructmat(:,chancmbindx) & selvec;
  chancmbneighbselmat(:,chancmbindx) = chancmbneighbselmat(:,chancmbindx) & selvec;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert Nx2 cell array with channel combinations into Nx1 cell array with labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = cmb2label(labelcmb);
label = {};
for i=1:size(labelcmb)
  label{i,1} = [labelcmb{i,1} '&' labelcmb{i,2}];
end


