function [params, s_new, avg, cnt] = denoise_avg2(params, s, state)
% DSS denoising function: Quasiperiodic averaging respecting 
%   trial boundaries in the original trials and with support for cell-arrays
%
%   [params, s_new, avg] = denoise_avg2(params, s, state)
%     params                Function specific modifiable parameters
%     params.tr             Trigger indices
%     params.tr_begin       How many samples before the trigger the ON state
%                           begins, if vectorial each element belongs to a 
%                           trigger
%     params.tr_end         How many samples after the trigger the ON state
%                           ends, if vectorial each element belongs to a trigger
%     params.tr_nsmp        Information pertaining to the number of samples of 
%                           each of the input trials (in the concatenated data-matrix)
%                           to avoid selection across trial boundaries
%     state                 DSS algorithm state
%     s                     Source signal estimate, matrix of row vector
%                           signals 
%     s_new                 Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: denoise_avg.m,v 1.19 2005/12/02 12:23:18 jaakkos Exp $

if nargin<3
  state = struct([]);
end

if ~isstruct(state)
    params.name        = 'Quasiperiodic averaging with known triggers';
    params.description = '';
    params.param       = {'fs','tr','tr_inds','begin','end'};
    params.param_value = {[], [], [], [], []};
    params.param_type  = {'scalar','vector','vector','scalar','scalar'};
    params.param_desc  = {'sampling frequency','trigger indices','used trigger indices','beginning of the ON mask','end of the ON mask'};
    params.approach    = {'pca','defl','symm'};
    params.alpha       = {};
    params.beta        = {'beta_global'};
    return;
end

if isfield(params, 'mask')
  mask = params.mask;
else
  mask = {};
end

% not available, using as direct indices
if isfield(params, 'tr_begin')
  tr_begin = params.tr_begin(:);
  tr_end   = params.tr_end(:);
end

% specify the indices for the triggering
if isfield(params,'tr_inds')
  tr_inds = params.tr_inds(:);
elseif ~isfield(params,'tr_inds') && isfield(params,'tr')
  tr_inds = params.tr(:);
else
  tr_inds = [];
end

% specify the behaviour of the demeaning, default is use the whole segment
demeanflag = 1;
if isfield(params, 'demean')
  if params.demean==0
    demeanflag = 0;
  elseif ischar(params.demean) && strcmp(params.demean, 'no')
    demeanflag = 0;
  elseif ischar(params.demean) && strcmp(params.demean, 'prezero')
    demeanflag = 2;
  elseif ischar(params.demean) && strcmp(params.demean, 'trialprezero')
    demeanflag = 3;
    if ~isfield(params, 'time'), error('demean with trialprezero only works with time axes in the params structure'); end
  end
end

% can be switched to false if the function is not used iteratively (for
% dss)
computenew = true;
if isfield(params, 'computenew') && params.computenew==0
  computenew = false;
end

if ~iscell(s)
  % input data matrix can be a concatenation of discontinuous trials, observe the edges
  if isfield(params, 'tr_nsmp') && numel(tr_begin)==1
    % update the tr_begin and tr_end scalars, to vectors
    
    % some quick sanity checks and conversions
    if sum(params.tr_nsmp(:))~=size(s,2)
      error('total number of samples according to params.tr_nsmp does not correspond with the size of the data');
    end
    if ~iscell(tr_inds)
      tr_inds = {tr_inds};
    end
    if ~iscell(tr_begin)
      tr_begin = {tr_begin};
    end
    if ~iscell(tr_end)
      tr_end = {tr_end};
    end
    
    [tr_inds, tr_begin, tr_end] = peaks2continuous(tr_inds, params.tr_nsmp, tr_begin, tr_end); 
    params.tr_inds  = tr_inds;
    params.tr_begin = tr_begin;
    params.tr_end   = tr_end;
  end
 
  if length(params.tr)==length(state.X)
    % trigger not as indices but as signal
    error('trigger indices must be extracted from the trigger signal first (automatic extraction not implemented yet)');
  end

  pst   = max(tr_end  - tr_inds);
  pre   = max(tr_inds - tr_begin);
  avg   = zeros(size(s,1),pre+pst+1);
  cnt   = zeros(size(avg));
  n     = tr_end-tr_begin+1;  

  % calculating the average
  for i = 1:length(tr_inds)
    begsmp = pre    - (tr_inds(i) - tr_begin(i)) + 1;
    endsmp = begsmp + (tr_end(i)  - tr_begin(i));
    tmp    = s(:,tr_begin(i):tr_end(i));
    if demeanflag,
      tmp    = tmp - sum(tmp,2)*(ones(1,size(tmp,2))./n(i));
    end
    avg(:,begsmp:endsmp) = avg(:,begsmp:endsmp) + tmp;
    cnt(:,begsmp:endsmp) = cnt(:,begsmp:endsmp) + 1;
  end
  avg = avg./cnt;
  
  % reconstructing the signals, only when needed
  if computenew,
    s_new = zeros(size(s));
    for i = 1:length(tr_inds)
      begsmp = pre    - (tr_inds(i) - tr_begin(i)) + 1;
      endsmp = begsmp + (tr_end(i)  - tr_begin(i));
      s_new(:,tr_begin(i):tr_end(i)) = s_new(:,tr_begin(i):tr_end(i)) + avg(:,begsmp:endsmp);
    end
  else
    s_new = [];
  end

elseif iscell(s)
  % ensure that pre and pst exists
  if ~isfield(params, 'pre')  params.pre = []; end
  if ~isfield(params, 'pst'), params.pst = []; end
  
  if isempty(tr_inds) && isfield(params, 'time')
    % lock to time point 0
    tr_inds = cell(1,numel(params.time));
    for k = 1:numel(params.time)
      tr_inds{k} = nearest(params.time{k},0);
    end
    params.tr_inds = tr_inds;
  end
    
  % input data cells are assumed to be discontinuous in time, observe the edges automatically
  if ~iscell(params.pre) && (numel(params.pre)==1 && numel(params.pst)==1)
    % scalar valued number of pre and post samples
    pre = tr_inds;
    pst = tr_inds;
    for k = 1:numel(tr_inds)
      pre{k}(:) = params.pre;
      pst{k}(:) = params.pst;
    end
    params.pre = pre;
    params.pst = pst;
  elseif ~iscell(params.pre) && isempty(params.pre) && isfield(params, 'time')
    pre = tr_inds;
    pst = tr_inds;
    for k = 1:numel(tr_inds)
      pre{k} = tr_inds{k}-1;
      pst{k} = numel(params.time{k})-tr_inds{k};
    end
    params.pre = pre;
    params.pst = pst;
  elseif iscell(params.pre)
    pre = params.pre;
    pst = params.pst;
  end 

  %maxpre = -inf;
  %maxpst = -inf;
  %for i = 1:length(pre)
  %  if ~isempty(pre{i}), maxpre = max([maxpre; pre{i}(:)]); end
  %  if ~isempty(pst{i}), maxpst = max([maxpst; pst{i}(:)]); end
  %end
  try
    maxpre = max(cat(2,pre{:}));
  catch
    maxpre = max(cat(1,pre{:}));
  end
  try
    maxpst = max(cat(2,pst{:}));
  catch
    maxpst = max(cat(1,pst{:}));
  end
  avg   = zeros(size(s{1},1),maxpre+maxpst+1);
  cnt   = zeros(1,           maxpre+maxpst+1);
  tim   = -maxpre:maxpst;
  
  % loop across the individual cells that contain useable data
  usetrials = find(~cellfun('isempty',tr_inds));
  for i = usetrials(:)'
    %if ~isempty(tr_inds{i})
      
      % create variables local to the i-loop
      dat   = s{i};
      [N,M] = size(dat);
      triggers = tr_inds{i};
      presmp   = pre{i};
      pstsmp   = pst{i};
      
      if ~isempty(mask)
        nonzero = mask{i}>0;
      else
        nonzero = ones(1,size(dat,2))>0;
      end
      begsmp = max(1, triggers-presmp);
      endsmp = min(M, triggers+pstsmp);
      
      if demeanflag==3,
        if ~exist('bslcnt', 'var')
          bslcnt = 0;
          bsl    = zeros(size(dat,1),1);
        end
        % keep track of the baseline
        bslsmp = find(nonzero(1:nearest(params.time{i},0)));
        bslcnt = bslcnt+numel(bslsmp);
        bsl    = bsl+sum(dat(:,bslsmp),2);
      end
      
      % loop across the triggers for the given cell
      for k = 1:numel(triggers)
        tmp        = dat(:,begsmp(k):endsmp(k));
        tmpnonzero = nonzero(begsmp(k):endsmp(k));
        %dat(:,~tmpnonzero) = nan;
        switch demeanflag
          case 1
            sel    = tmpnonzero;
            n      = sum(sel);
            if n==0, continue; end
          case 2
            sel    = tmpnonzero(1:min(presmp(k),size(tmp,2)));
            n      = sum(sel);
            if n==0, continue; end
        end
        
        if demeanflag>0 && demeanflag<3,
          %tmp = tmp - (sum(tmp(:,sel),2)/n)*ones(1,size(tmp,2));
          tmp = bsxfun(@minus, tmp, sum(tmp(:,sel),2)/n);
        end
        
        % get the indices to the samples in the average that have data
        % for the current snippet of data
        tmpindx = tim>=begsmp(k)-triggers(k) & tim<=endsmp(k)-triggers(k);
        
        % do a running sum
        avg(:,tmpindx) = avg(:,tmpindx) + tmp;
        cnt(1,tmpindx) = cnt(1,tmpindx) + tmpnonzero;
      end % for k
      
    %end % if 
  end % for i
  
  % normalize
  avg = avg./cnt(ones(N,1),:);
  
  if demeanflag==3,
    avg = avg - repmat(bsl./bslcnt, [1 size(avg,2)]);
  end
  
  if computenew
    s_new = cell(size(s));
    for i = 1:length(s)
      s_new{i} = zeros(size(s{i}));
      if ~isempty(tr_inds{i})
        for k = 1:length(tr_inds{i})
          begsmp = max(1,tr_inds{i}(k)-pre{i}(k));
          endsmp = min(size(s_new{i},2),tr_inds{i}(k)+pst{i}(k));
          tmpindx = tim>=begsmp-tr_inds{i}(k) & tim<=endsmp-tr_inds{i}(k);
          s_new{i}(:, begsmp:endsmp) = s_new{i}(:, begsmp:endsmp) + avg(:, tmpindx);
        end
      end
    end
  else
    s_new = cell(size(s));
  end
  
end

function [p, begsmp, endsmp] = peaks2continuous(p, nsmp, presmp, postsmp)

% convert cell-array with peak indices into a vector, where the indices
% are adjusted according to the nsmp per trial, and the begin and endsamples
% respect the trial boundary.
% 
% Use as
%   [p, begsmp, endsmp] = peaks2continuous(p, sampleinfo, presmp, postsmp) 

ix = zeros(0,1);
iy = zeros(0,1);
for k = 1:numel(p)
  ix = [ix; k.*ones(numel(p{k}),1)];
  iy = [iy; p{k}(:)];
end

if numel(presmp)==1 && numel(presmp{1})==1
  presmp = repmat(presmp{1}(1),[numel(ix) 1]);
else
  presmp = cell2mat(presmp)';
end
if numel(postsmp)==1 && numel(postsmp{1})==1
  postsmp = repmat(postsmp{1}(1),[numel(ix) 1]);
  postsmp = cell2mat(postsmp)';
end

csmp   = cumsum([0 nsmp(:)']);
begsmp = zeros(numel(ix),1);
endsmp = zeros(numel(ix),1);
for k = 1:numel(ix)
  begsmp(k,1) = csmp(ix(k)) + max(iy(k)-presmp(k),  1);
  endsmp(k,1) = csmp(ix(k)) + min(iy(k)+postsmp(k), nsmp(ix(k)));
  iy(k)       = csmp(ix(k)) + iy(k);
end
p = iy;

