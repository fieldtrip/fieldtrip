function [metric,all] = evaluate(post,tcls,varargin)
%EVALUATOR evaluation criterion for classifiers/regressors
%
%   [metric,all] = evaluate(post,tcls,varargin)
%
%   metric returns an average metric and all the metrics of all inputs if
%   the input is a cell array
%
%   parameter 'metric' determines the evaluation criterion:
%
%   classifiers:
%       'accuracy' : classification accuracy
%       {'meanrate';'maxrate';'minrate';'medianrate';'sumrate'} :
%         functions that can be applied to the diagonal of a confusion matrix
%       'tabcounts' : contingency table
%       'cfmatrix' : confusion matrix
%       'kappa' : kappa statistic
%       'chi2counts' : chi-square statistic for contingency table
%       'negloglik' : negative log likelihood of posteriors
%       'loglik' : log likelihood of posteriors
%       'dprime' : the difference between z-transformed hits and false-alarms
%       'auc' : area under the receiver-operating characteristic curve
%       'raw' : raw results; posterior and true classes
%       'assignments' : class assignments
%       'mi' : mutual information under assumption of uniform class priors
%               or using the probabilities defined as 'prior'
%       'itr' : information transfer rate (bits per minute); requires 'duration' of
%               a classification be defined in seconds.
%       'tp' ('hits'): true positives
%       'fp' : false positives
%       'tn' : true negatives
%       'fn' ('misses'): false negatives
%       'precision' : 1 - FP/(TP + FP)
%       'recall' ('sensitivity','hitrate'): TP/(TP + FN)
%       'missrate' (a.k.a. false detection rate) : FP/(TP + TP)
%       'specificity' : TN/(TN + FP)
%       'f1' : 2*TP/(2*TP + FN + FP)
%       'hfdiff' : hitrate - missrate (or sensitivity + precision - 1)
%
%   regressors:
%       'rss' : residual sum of squares
%       'circrss' : circular residual sum of squares
%
%   Copyright (c) 2008, Marcel van Gerven, Christian Hesse
%
%   $Log: evaluate.m,v $
%

  options = varargin2struct(varargin);
  
  if ~isfield(options,'metric'), options.metric = 'accuracy'; end

  % compute metric for all cell elements
  if iscell(post)

    warning off all
    all = cell(1,length(post));
    metric = cell(1,length(post));

    for c=1:length(post)
      [metric{c},all{c}] = evaluate(post{c},tcls{c},varargin{:});
    end
    warning on all

    % concatenate all results and compute metric
    apost = [];
    atcls = [];
    for c=1:length(post)
      apost = cat(1,apost,post{c});
      atcls = cat(1,atcls,tcls{c});
    end

    metric = evaluate(apost,atcls,varargin{:});

  else

    metric = compute_metric(post,tcls(:,1),options);
    all = metric;
  end
end

function met = compute_metric(post,tcls,cfg)

% precompute confusion matrix as it is often used
[cobs,cexp] = contingency(post,tcls);

nclasses = size(post,2);

met = [];
switch lower(cfg.metric)

  case 'accuracy'

    % compute the total correct classification rate which collapses over all
    % classes and is insensitive to the relative frequency of occurrence of
    % examples from each class
    met = sum(diag(cobs))./size(tcls,1);

    % NOTE: if the number of samples per class is the same then this metric
    % is the same as 'meanrate' in the next option, otherwise the measure is
    % dominated by the class(es) with the most samples

  case {'meanrate';'maxrate';'minrate';'medianrate';'sumrate'}

    % get the name of the function to apply to the diagonal
    tmp = lower(cfg.metric);
    tmp = tmp(1:end-length('rate')); % <- this is just to be explicit
    cmd = ['met = ',tmp,'(diag(contingency(post,tcls)));'];
    eval(cmd);

    % met should now be appropriately assigned

  case 'tabcounts'

    % return contingency table
    met = cobs;

  case 'cfmatrix'

    % return the table of classification rates, i.e.,
    % the so-called confusion matrix

    met = cobs ./ repmat(sum(cobs,2),[1 size(cobs,2)]);
    met(isnan(met)) = 0;

  case 'kappa'

    % return Cohen's kappa coefficient

    % number of trials
    N = sum(cobs(:));

    % compute classification accuracy
    p0 = sum(diag(cobs))/N;

    % compute chance agreement
    pe = 0;
    for i=1:size(cobs,1)
      pe = pe + sum(cobs(i,:))*sum(cobs(:,i));
    end
    pe = pe / N^2;

    % compute kappa metric
    met = (p0 - pe)/(1 - pe);


  case 'chi2counts'

    % compute the chi-square statistic for the contingency table
    met = sum( ((cobs(:) - cexp(:)).^2) ./ cexp(:) );

  case {'negloglik'; 'loglik'}

    % compute log-likelihood of the real class labels
    tpost = post';
    met = sum(log(tpost((0:size(tpost,1):((size(tpost,2)-1)*size(tpost,1))) + tcls')));

    if strcmpi(cfg.metric,'negloglik')
      met = - met;
    end

  case 'lik'

    % compute likelihood of the real class labels
    tpost = post';
    met = prod(tpost((0:size(tpost,1):((size(tpost,2)-1)*size(tpost,1))) + tcls'));

  case 'dprime'

    % note that this is only defined for 2 classes
    if nclasses > 2
      error(['option ''',cfg.metric,''' only valid for 2 classes']);
    end

    % compute the difference between z-transformed hits and false-alarms (this
    % can lead to +/- Inf if one of the rates is 0 or 1)
    met = norminv(cobs(1,1))-norminv(cobs(1,2));

    % NOTE: this requires the statistics toolbox

  case 'auc'

    % area under the receiver-operating characteristic curve

    % note that this is only defined for 2 classes
    if nclasses > 2
      error(['option ''',cfg.metric,''' only valid for 2 classes']);
    end

    ppos = post(tcls == 1,1);
    pneg = post(tcls == 2,1);

    met = 0;
    for i=1:length(ppos)
      for j=1:length(pneg)

        if ppos(i) > pneg(j)
          met = met + 1;
        elseif ppos(i) == pneg(j)
          met = met + 0.5;
        end
      end
    end
    met = met / (length(ppos) * length(pneg));

  case 'raw'

    % return raw results
    met = [post tcls];

  case 'assignments'

    % return for each case if the correct class is assigned; needed e.g. for binomial test
    met = (pcls == tcls);

  case {'mi' 'itr'}

    % normalize prior if specified
    if isfield(cfg,'prior')
      cfg.prior = cfg.prior./sum(cfg.prior);
    end

    % confusion matrix
    cfm = cobs ./ repmat(sum(cobs,2),[1 size(cobs,2)]);
    cfm(isnan(cfm)) = 0;

    % compute mutual information
    met = 0;
    M = size(cobs,1);
    for i=1:M

      if isfield(cfg,'prior')
        % use prespecified prior
        prior = cfg.prior(i);
      else
        % assuming uniform class priors
        prior = 1/M;
      end

      for j=1:M
        if cobs(i,j)
          pyx = cfm(i,j)/sum(cfm(i,:));
          met = met + prior * pyx * log2(pyx);
        end
      end

    end

    for j=1:M

      py = sum(cfm(:,j))/M;
      if py
        met = met - py*log2(py);
      end
    end
    
    met = max(met,0);

    if strcmpi(cfg.metric,'itr')

      if ~isfield(cfg,'duration'), error('duration (s) of a classification needs to be specified for ITR'); end

      % duration of a classification in seconds
      % computes information transfer rate in bits per minute

      met = met * 60/cfg.duration;

    end

  case {'tp','hits','precision','recall','sensitivity','f1','hfdiff'}

    if nclasses == 2 % for binary, report true positives w.r.t. class 1
      tp = cobs(1,1);
    else % for multiclass, report true positives for each class
      tp = diag(cobs);
    end
    met = tp;

  case {'fp','precision','specificity','f1','hfdiff'}

    if nclasses == 2 % for binary, report false positives w.r.t. class 1
      fp = cobs(2,1);
    else % for multiclass, report false positives for each class
      fp = sum(cobs,1)' - diag(cobs);
    end
    met = fp;

  case {'tn','specificity'}

    if nclasses == 2 % for binary, report true negatives w.r.t. class 1
      met = cobs(2,2);
    else % for multiclass, report true negatives for each class
      met = zeros(1,nclasses);
      for j=1:nclasses
        met = sum(cobs(:)) - sum(cobs(j,:)) - sum(cobs(:,j));
      end
    end

  case {'fn','misses','recall','sensitivity','f1','hfdiff'}

    if nclasses == 2 % for binary, report false negatives w.r.t. class 1
      met = cobs(1,2);
    else % for multiclass, report false negatives for each class
      met = sum(cobs,2) - diag(cobs);
    end

  case 'precision'

    met = 1 - (fp/(tp+fp));

  case {'recall','sensitivity','hitrate'}

    met = tp/(tp+fn);

  case 'specificity'

    met = tn/(tn+fp);

  case 'missrate' % false detection rate

    met = fp/(tp + fp);

  case 'f1'

    met = 2*tp/(2*tp + fn + fp);

  case 'hfdiff'

    met = (tp/(tp+fn)) - (fp/(tp+fp));

    % regressor metrics; assume the first column of post is the predicted
    % value

  case 'rss' % residual sum of squares

    met = sum((tcls - post(:,1)).^2);

  case 'circrss' % circular residual sum of squares; assumes data is in the range -pi .. pi

    met = mod(mod(tcls,2*pi) - mod(post(:,1),2*pi),2*pi);
    met = min(met,2*pi - met);
    met = sum(met.^2);

  otherwise
    error(['unsupported option ',cfg.metric]);

end

end

function [cobs,cexp] = contingency(post,tcls)
% compute contingency table with true class as rows and predicted class as
% columns

  % number of classes
  numcls = size(post,2);

  % predicted classes with the maximum posterior
  % probability
  [temp,pcls] = max(post,[],2);

  cobs = zeros(numcls);
  cexp = zeros(numcls);
  for i=1:numcls % true class is in rows
    ii = (tcls==i);
    for j=1:numcls % assigned class is in columns
      jj = (pcls==j);
      cobs(i,j) = sum(ii(:) & jj(:));
      cexp(i,j) = sum(ii(:))/numcls;

    end

  end

end
