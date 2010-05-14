classdef classifier < predictor
%CLASSIFIER abstract classifier method class
%
% A classifier takes a variable number of arguments upon construction. 
% During operation, the classifier takes data and
% produces posterior probabilities of class labels as an N x C matrix for N
% examples and C classes.
% 
% Subclasses should implement the train and test functions and optionally
% the getmodel function which reshapes parameters to something
% interpretable
%
% SEE ALSO:
% doc classifiers
%
% Copyright (c) 2008, Marcel van Gerven


  methods
  
    function obj = classifier(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      
    end
    
    function clf = predict(obj,data)
      % convert posterior dataset into classifications
      
      [tmp,clf] = max(obj.test(data),[],2);
    end
    
  end
  
  methods(Static)
    
    function metric = evaluate(post,tcls,varargin)
      %EVALUATE evaluation metrics for classifiers
      %
      %   metric = evaluate(post,tcls,varargin)
      %
      %   parameter 'metric' determines the evaluation criterion:
      %
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
      %               a classification to be defined in seconds.
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
      
      options = struct(varargin{:});
      
      if ~isfield(options,'metric'), options.metric = 'accuracy'; end
      
      metric = compute_metric(post,tcls(:,1),options);
      
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
      
    end
    
    function p = significance(post,design,varargin)
      % Compute significance level that our result is different from
      % another classifier.
      %
      % tests:
      %
      % 'binomial_test' (default)
      % Comparison is based on a one-sided binomial test (McNemar)
      % Bonferroni correction is applied for multiple datasets
      %
      
      options = struct(varargin{:});
      
      p = cell(1,length(post));
      for c=1:length(post)
        
        if ~isfield(options,'test') || strcmp(options.test,'binomial_test')
          % one-sided binomial test with automatic bonferroni correction
          
          % compute class with highest prior probability
          nclasses = size(post{c},2);
          priors = zeros(1,nclasses);
          for k=1:nclasses
            priors(k) = priors(k) + sum(design{c}(:,1)==k);
          end
          [mxx,clss] = max(priors);
          
          rndpost = zeros(size(post{c}));
          rndpost(:,clss) = 1;
          
          [r,p{c},level] = classifier.binomial_test(post{c},design{c},rndpost,design{c},'twosided',false,'bonferroni',length(post));
          
          if r
            fprintf('dataset %d: null hypothesis rejected (%g<%g);\nsignificant difference from ',c,p{c},level);
          else
            fprintf('dataset %d: null hypothesis not rejected (%g>%g);\nno significant difference from ',c,p{c},level);
          end
          fprintf('majority classification (class %d with prior of %g)\n',clss,mxx/sum(priors));
          
        end
      end
      
      if length(p) == 1
        p = p{1};
      end
      
    end
    
    function [reject,pvalue,level] = binomial_test(cpost1, cdesign1, cpost2, cdesign2, varargin)
      % BINOMIAL_TEST makes a significance test whether two algorithms perform
      % the same or differently
      %
      % [reject,pvalue,level] = significance_test(cpost1, design1, cpost2, design2, varargin)
      %
      % where typically design1 == design2
      %
      % reject = reject the null hypothesis? I.e., is there a difference?
      % pvalue = used pvalue
      % level = used (Bonferroni corrected) significance level
      %
      % 'level' = 0.05; signifance level for rejecting the null hypothesis
      % 'test' = 'mcnemar'; type of significance test
      % 'twosided' = true; one versus two-sided test
      % 'bonferroni' = 1; lists number of tests;
      %   if larger than 1 we do bonferroni correction for multiple tests
      % 'idx'; gives the column index of the design matrix that represents
      %    example indices. Used when design1 and design2 are not the same
      %
      % if inputs are cell arrays then all experiments are
      % concatenated (e.g., folds)
      
      % initialization
      
      cfg = struct(varargin{:});
      
      if ~isfield(cfg,'test'), cfg.test = 'mcnemar'; end
      if ~isfield(cfg,'bonferroni'), cfg.bonferroni = 1; end
      if ~isfield(cfg,'level'), level = 0.05; else level = cfg.level; end
      if ~isfield(cfg,'twosided'), cfg.twosided = true; end
      
      % change significance level
      if cfg.bonferroni > 1, level = 1 - power(1 - level,1 / cfg.bonferroni); end
      
      % process cell arrays
      
      if iscell(cpost1)
        post1 = [];
        for c=1:numel(cpost1)
          post1 = [post1; cpost1{c}];
        end
        cpost1 = post1;
      end
      
      if iscell(cpost2)
        post2 = [];
        for c=1:numel(cpost2)
          post2 = [post2; cpost2{c}];
        end
        cpost2 = post2;
      end
      
      % true classes
      
      if iscell(cdesign1)
        tcls = [];
        for c=1:numel(cdesign1)
          tcls = [tcls; cdesign1{c}];
        end
        ctcls1 = tcls;
      else
        ctcls1 = cdesign1;
      end
      
      if iscell(cdesign2)
        tcls = [];
        for c=1:numel(cdesign2)
          tcls = [tcls; cdesign2{c}];
        end
        ctcls2 = tcls;
      else
        ctcls2 = cdesign2;
      end
      
      % optionally shuffle rows of the design matrices
      if isfield(cfg,'idx')
        
        [a,b] = sort(ctcls1(:,cfg.idx));
        ctcls1 = ctcls1(b,:);
        cpost1 = cpost1(b,:);
        
        [a,b] = sort(ctcls2(:,cfg.idx));
        ctcls2 = ctcls2(b,:);
        cpost2 = cpost2(b,:);
        
        % now both designs should be comparable
      else
        
        % check if designs are the same if index is left unspecified
        if any(ctcls1(:,1) ~= ctcls2(:,1))
          error('designs do not match; please specify indices');
        end
      end
      
      % execute a one-sided test
      switch cfg.test
        
        case 'mcnemar'
          pvalue = mcnemar(cpost1,cpost2,ctcls1);
          
      end
      
      % make two-sided if necessary
      if cfg.twosided, level = level/2; end
      
      % make decision if we can reject the null hypothesis
      reject = (pvalue < level);
      
      function pvalue = mcnemar(cpost1, cpost2, ctcls2)
        % MCNEMAR tests whether the posterior probabilities cpost1 and cpost2
        % are significantly different.
        %
        % Note: this is an approximation to the computationally heavy binomial
        % test. Later, this file should become a general significance testing
        % function that also accepts cell arrays for cpost1 and cpost2
        
        % compute class labels
        [temp,pcls1] = max(cpost1,[],2);
        [temp,pcls2] = max(cpost2,[],2);
        
        % compute winners for each algorithm
        
        p1 = (pcls1 == ctcls2(:,1));
        p2 = (pcls2 == ctcls2(:,1));
        
        s = sum(p1 & ~p2); % successes w.r.t. 1
        f = sum(~p1 & p2); % failures w.r.t. 1
        
        % compute one sided p-values
        
        if (s+f)
          mcnemarstat = (abs(s - f) - 1).^2 / (s+f);
        else
          mcnemarstat = 0;
        end
        
        pvalue = 1 - chi2cdf(mcnemarstat,1);
        
      end
      
    end
    
  end
  
end
