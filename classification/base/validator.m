classdef validator
  %VALIDATOR abstract validator class
  %
  %   a validator takes a classification procedure
  %
  %   data is not randomized by default by a validator since it may have an inherent
  %   temporal ordering (i.e, for dynamic models)
  %
  %   Options:
  %   'init'      : initialize the random number generator with seed [1]
  %   'randomize' : reshuffle if true [false]
  %   'procedure' : the classification procedure
  %   'verbose'   : output comment if true [false]
  %   'balanced'  : balance data before training/testing a classifier (true,
  %                   false,'train', 'test'); [false]
  %   'compact'   : only retain posteriors and design (in case of memory issues)
  %   'mode'      : 'classification' or 'regression' []
  %   'transfer'  : transfer learning []
  %
  %   SEE ALSO:
  %   crossvalidator.m
  %   loocrossvalidator.m
  %
  %   Copyright (c) 2008, Marcel van Gerven
  %
  %   $Log: validator.m,v $
  %
  
  properties
    
    verbose = false;
    
    post;   % posteriors
    design  % associated class labels
    
    init        = 1;     % initialize RNG
    randomize   = false; % randomize trials
    balanced    = false; % balance classes
    compact     = false; % only retain necessary results
    
    % operating mode (classification/regression);
    % for clf, classes are evenly distributed over folds
    mode;
    
    % transfer learning
    transfer; 
    
    procedure;
    
  end
  
  methods
    function obj = validator(varargin)
      % a procedure is a clfproc or a cell array of clfprocs
      % the latter is useful for retaining all results in a
      % crossvalidation (although it takes up more space)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      
      if isempty(obj.procedure)
        error('procedure must be specified!');
      end
      
      if ~isa(obj.procedure,'clfproc')
        % try to create procedure if input is a cell array or a predictor
        if obj.verbose
          fprintf('creating procedure\n');
        end
        obj.procedure = clfproc(obj.procedure);
      end
      
      % determine mode of operation
      if isempty(obj.mode)
        
        m = obj.getpredictor();
        
        if isa(m,'classifier') || (isa(m,'optimizer') && isa(m.method,'classifier'))
          obj.mode = 'classification';
        elseif isa(m,'regressor') || (isa(m,'optimizer') && isa(m.method,'regressor'))
          obj.mode = 'regression';
        else
          obj.mode = nan;
        end
        
      end
      
      % determine transfer learning
      if isempty(obj.transfer)
        
        m = obj.getpredictor();
        
        if isa(m,'transfer_learner') || (isa(m,'optimizer') && isa(m.method,'transfer_learner'))
          obj.transfer = true;
        else
          obj.transfer = false;
        end
        
      end
      
      if obj.verbose
        fprintf('creating validator for clfproc %s\n',obj.procedure.name);
      end
      
    end
    
    function n = nclasses(obj)
      % return number of classes when known (called by statistics_crossvalidate)
      
      n = obj.getpredictor().nclasses;
            
    end        
        
    function m = getmodel(obj,label,dims)
      % try to return the classifier parameters as a model
      % wrt some class label and reshape into dims when specified
      
      if nargin < 2, label = 1; end      
      if nargin < 3, dims = []; end
            
      if ~iscell(obj.procedure)
        obj.procedure = {obj.procedure};
      end
      
      if isempty(obj.procedure{1})
        m = [];
        return;
      end
      
      fm = cellfun(@(x)(x.getmodel(label,dims)),obj.procedure,'UniformOutput', false);
      
      m = fm{1};
      for c=2:length(obj.procedure)
        
        if iscell(m)
          % if we return multiple models
          
          for j=1:length(m)
            m{j} = m{j} + fm{c}{j};
          end
          
        else
          m = m + fm{c};
        end
      end
      
      % take the mean of the parameters
      if iscell(m)
        for c=1:length(m)
          m{c} = m{c}./length(obj.procedure);
        end
      else
        m = m./length(obj.procedure);
      end
            
    end
    
    
    function [result,all] = evaluate(obj,varargin)
      % indirect call to evaluate to simplify the interface and handle
      % transfer learned data
      
      % check if we are dealing with transfer learning
      if obj.transfer
        
        [res all] = validator.eval(obj.post,obj.design,varargin{:});
        
        % generate the result for all subjects
        if iscell(obj.post) && iscell(obj.post{1})
          % n-fold result
          
          result = zeros(1,length(obj.post{1}));
          for k=1:length(obj.post{1})
            for j=1:length(obj.post)
              result(k) = result(k) + all{j}{k};
            end
          end
          result = result./length(obj.post);
          
        else
          % percentage
          
          result = zeros(1,length(obj.post));
          for k=1:length(obj.post)
            result(k) = all{k};
          end
          
        end
        
      else        

        [result all] = validator.eval(obj.post,obj.design,varargin{:});
        
      end
      
    end    
    
    function p = significance(obj)
      % Compute significance level that our result is different from
      % another classifier.
      %
      % Comparison is based on a one-sided binomial test (McNemar)
      % without Bonferroni correction
      %
      
      % check if we are dealing with transfer learning
     if obj.transfer
     
        tmppost = obj.post;
        tmpdesign = obj.design;
     
        % generate the result for all subjects
        if iscell(tmppost) && iscell(tmppost{1})
          % n-fold result
                     
          p = zeros(1,length(tmppost{1}));
          
          for k=1:length(tmppost{1})
            
            tpost = cell(1,length(tmppost));
            tdesign = cell(1,length(tmppost));
            for j=1:length(tmppost)
              tpost{j}   = tmppost{j}{k};
              tdesign{j} = tmpdesign{j}{k};
            end
            
            p(k) = obj.compute_significance(tpost,tdesign);
            
          end
          
          
        else
          % percentage                    
          
          p = zeros(1,length(tmppost));
          
          for k=1:length(tmppost)
            
            tpost   = tmppost{k};
            tdesign = tmpdesign{k};
            
            p(k) = obj.compute_significance(tpost,tdesign);
          end
          
        end
                  
      else
        
        p = obj.compute_significance(obj.post,obj.design);
        
      end
      
    end
  end
  
  methods(Access = protected)
    
    function p = compute_significance(obj,tpost,tdesign)
              
      if iscell(tpost)
        
        % compute class with highest prior probability
        nclasses = size(tpost{1},2);
        priors = zeros(1,nclasses);
        for c=1:length(tpost)
          for k=1:nclasses
            priors(k) = priors(k) + sum(tdesign{c}(:,1)==k);
          end
        end
        [mxx,class] = max(priors);
        
        rndpost = cell(1,length(tpost));
        for c=1:length(tpost)
          rndpost{c} = zeros(size(tpost{c}));
          rndpost{c}(:,class) = 1;
        end
      else
        
        % compute class with highest prior probability
        nclasses = size(tpost,2);
        priors = zeros(1,nclasses);
        for k=1:nclasses
          priors(k) = sum(tdesign(:,1)==k);
        end
        [mxx,class] = max(priors);
        
        rndpost = zeros(size(tpost));
        rndpost(:,class) = 1;
      end
      
      if obj.verbose
        fprintf('performing one-sided biniomial test (p=0.05)\n');
      end
      [r,p,level] = validator.significance_test(tpost,tdesign,rndpost,tdesign,'twosided',false);
      
      if obj.verbose
        
        if r
          fprintf('null hypothesis rejected (%g<%g);\nsignificant difference from ',p,level);
        else
          fprintf('null hypothesis not rejected (%g>%g);\nno significant difference from ',p,level);
        end
        
        if (nargin == 2 && random)
          fprintf('random classification.\n');
        else
          fprintf('majority classification (class %d with prior of %g)\n',class,mxx/sum(priors));
        end
      end
    end      
    
    function [newdata,newdesign] = balance(obj,data,design)
      % balance data; make sure classes are evenly represented by
      % resampling with replacement
      
      if iscell(data)
        
        newdata = cell(1,length(data));
        newdesign = cell(1,length(design));
        for c=1:length(data)
            [newdata{c},newdesign{c}] = obj.balance(data{c},design{c});
        end
        
      else
        
        nclasses = max(design(:,1));
        maxsmp = 0;
        for j=1:nclasses
          summed = sum(design(:,1) == j);
          if ~summed
            fprintf('not all classes represented in data; refusing to balance\n');
            newdata = data;
            newdesign = design;
            return;
          end
          maxsmp = max(maxsmp,summed);
        end
        newdata = zeros(nclasses*maxsmp,size(data,2));
        for j=1:nclasses
          cdata = data(design(:,1) == j,:);
          if  size(cdata,1) ~= maxsmp
            
            % sample with replacement
            newdata(((j-1)*maxsmp+1):(j*maxsmp),:) = cdata(ceil(size(cdata,1)*rand(maxsmp,1)),:);
          else
            newdata(((j-1)*maxsmp+1):(j*maxsmp),:) = cdata;
          end
        end
        newdesign=ones(nclasses*maxsmp,1);
        for j=2:nclasses
          newdesign(((j-1)*maxsmp+1):(j*maxsmp)) = j;
        end
        
        % shuffle data
        [newdata,newdesign] = obj.shuffle(newdata,newdesign);
        
      end
    end
    
    function data = collapse(obj,data)
      % check dimensions of data
      
      if obj.verbose, fprintf('collapsing data\n'); end
      
      if iscell(data)
        for c=1:length(data)
          if ndims(data{c}) > 2
            
            sz = size(data{c});
            data{c} = reshape(data{c},[sz(1) prod(sz(2:end))]);
          end
        end
      else
        if ndims(data) > 2
          
          sz = size(data);
          data = reshape(data,[sz(1) prod(sz(2:end))]);
        end
      end
    end
    
    function [data,design] = shuffle(obj,data,design)
      % shuffle (randomize) examples
      %
      %   Copyright (c) 2008, Marcel van Gerven
      %
      %   $Log: shuffle.m,v $
      %
      
      if obj.verbose, fprintf('shuffling data\n'); end
      
      if iscell(data)
        
        % try to keep the same permutation when possible
        prm = randperm(size(data{1},1))';
        
        if ~iscell(design), design = design(prm,:); end
        
        for c=1:length(data)
          
          sz = size(data{c});
          
          if size(prm,1) ~= sz(1)
            prm = randperm(sz(1));
          end
          
          % first randomize the ordering of the data
          data{c} = reshape(data{c}(prm,:),sz);
          
          if iscell(design)
            design{c} = design{c}(prm,:);
          end
        end
      else
        sz = size(data);
        
        prm = randperm(sz(1));
        data = reshape(data(prm,:),sz);
        design = design(prm,:);
      end
    end
    
    function p = getpredictor(obj)
   
      p = [];
      if iscell(obj.procedure)
        if ~isempty(obj.procedure{1})
          p = obj.procedure{1}.clfmethods{end};
        end
      else
        if ~isempty(obj.procedure)
          p = obj.procedure.clfmethods{end};
        end
      end
      
      if obj.verbose && isempty(p)
        fprintf('could not determine classifier type\n');
      end
      
    end
    
    function [ldat,ldes,udat,udes] = get_labeled_unlabeled(obj,data,design)
      % return labeled and unlabeled parts; nan indicates missing label

      if iscell(data)
        
        udat = cell(1,length(data));
        ldat = cell(1,length(data));
        udes = cell(1,length(data));
        ldes = cell(1,length(data));
        for c=1:length(data)
          
          udes{c} = isnan(design{c});
          ldes{c} = ~isnan(design{c});
          
          udat{c} = data{c}(udes{c},:);
          ldat{c} = data{c}(ldes{c},:);
          
          udes{c} = nan(nnz(udes{c}),1);
          ldes{c} = design{c}(ldes{c});
          
        end
        
      else
        
        udes = isnan(design);
        ldes = ~isnan(design);
        
        udat = data(udes,:);
        ldat = data(ldes,:);
        
        udes = nan(nnz(udes),1);
        ldes = design(ldes);        
        
      end
        
    end
    
  end
  
  methods(Static)
    
    function [metric,all] = eval(post,tcls,varargin)
      %EVAL evaluation criterion for classifiers/regressors
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
          [metric{c},all{c}] = validator.eval(post{c},tcls{c},varargin{:});
        end
        warning on all
        
        % concatenate all results and compute metric
        apost = [];
        atcls = [];
        for c=1:length(post)
          apost = cat(1,apost,post{c});
          atcls = cat(1,atcls,tcls{c});
        end
        
        metric = validator.eval(apost,atcls,varargin{:});
        
      else
        
        metric = compute_metric(post,tcls(:,1),options);
        all = metric;
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
      
    end
    
    
    function [reject,pvalue,level] = significance_test(cpost1, design1, cpost2, design2, varargin)
      % SIGNIFICANCE_TEST makes a significance test whether two algorithms perform
      % the same or differently
      %
      % [reject,pvalue,level] = significance(cpost1, design1, cpost2, design2, varargin)
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
      %
      % Copyright (C) 2008, Marcel van Gerven
      % F.C. Donders Centre for Cognitive Neuroscience, Nijmegen, NL
      %
      
      % initialization
      
      cfg = varargin2struct(varargin);
      
      if ~isfield(cfg,'test'), cfg.test = 'mcnemar'; end
      if ~isfield(cfg,'bonferroni'), cfg.bonferroni = 1; end
      if ~isfield(cfg,'level'), level = 0.05; else level = cfg.level; end
      if ~isfield(cfg,'twosided'), cfg.twosided = true; end
      
      % change significance level
      if cfg.bonferroni > 1, level = 1 - power(1 - level,1 / cfg.bonferroni); end
      
      % process cell arrays
      
      if iscell(cpost1)
        post1 = [];
        for c=1:length(cpost1)
          post1 = [post1; cpost1{c}];
        end
        cpost1 = post1;
      end
      
      if iscell(cpost2)
        post2 = [];
        for c=1:length(cpost2)
          post2 = [post2; cpost2{c}];
        end
        cpost2 = post2;
      end
      
      % true classes
      
      if iscell(design1)
        tcls = [];
        for c=1:length(design1)
          tcls = [tcls; design1{c}];
        end
        ctcls1 = tcls;
      else
        ctcls1 = design1;
      end
      
      if iscell(design2)
        tcls = [];
        for c=1:length(design2)
          tcls = [tcls; design2{c}];
        end
        ctcls2 = tcls;
      else
        ctcls2 = design2;
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
