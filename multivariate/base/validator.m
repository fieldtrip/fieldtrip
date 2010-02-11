classdef validator
    
  properties
    
    procedure;
    
    post;   % posteriors
    design  % associated class labels
    
    init        = 5;    % initialize RNG
    balanced    = false; % balance classes
    compact     = false; % only retain necessary results
    verbose     = false; % verbose output
    
    
  end
  
  methods
    function obj = validator(varargin)
      % a procedure is a mva
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
      
      if isempty(obj.procedure)
        error('procedure must be specified!');
      end
      
      if ~isa(obj.procedure,'mva')
        % try to create procedure if input is a cell array or a predictor
        obj.procedure = mva(obj.procedure);
      end
      
      if obj.verbose
        fprintf('creating validator for mva %s\n',obj.procedure.name);
      end
      
      % initialize random number generator
      if ~isempty(obj.init)
        if obj.verbose
          fprintf('initializing random number generator with seed %d\n',obj.init)
        end
        RandStream.setDefaultStream(RandStream('mt19937ar','seed',obj.init));
      end
      
    end
    
%     function n = nclasses(obj)
%       % return number of classes when known (called by statistics_crossvalidate)
%       
%       n = obj.getpredictor().nclasses;
%       
%     end
    
    function [m,desc] = getmodel(obj)
      % try to return the classifier parameters as a model
          
      if ~iscell(obj.procedure)
        obj.procedure = {obj.procedure};
      end
      
      if isempty(obj.procedure{1})
        
        m = {};
        desc = {};
        return;
      
      end
      
      fm = cellfun(@(x)(x.getmodel()),obj.procedure,'UniformOutput', false);
      
      [tmp,desc] = obj.procedure{1}.getmodel();
      
      m = fm{1};
      for c=2:length(obj.procedure)
        
        for j=1:length(m)
          m{j} = m{j} + fm{c}{j};
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
    
    
    function result = evaluate(obj,varargin)
      % indirect call to eval to simplify the interface
      
      % create the concatenation of all folds for each of the datasets
     
      post = obj.post;
      design = obj.design;
      
      tpost = cell(1,size(post,2));
      tdesign = cell(1,size(design,2));
      for c=1:length(tpost)
        tpost{c} = cat(1,post{:,c});
        tdesign{c} = cat(1,design{:,c});
      end
      
      result = cell(1,length(tpost));
      
      for c=1:length(tpost)
        result{c} = obj.getpredictor().evaluate(tpost{c},tdesign{c},varargin{:});   
      end
      
      if length(result) == 1
        result = result{1};
      end
      
    end
    
    function p = significance(obj,varargin)
      % Compute significance level that our result is different from
      % another classifier.
      %
      % tests:
      %
      % 'binomial_test' (default)
      % Comparison is based on a one-sided binomial test (McNemar)
      % Bonferroni correction is applied for multiple datasets
      %
      options = varargin2struct(varargin);
      
      post = obj.post;
      design = obj.design;
      
      % create the concatenation of all folds for each of the datasets
      tpost = cell(1,size(post,2));
      tdesign = cell(1,size(design,2));
      for c=1:length(tpost)
        tpost{c} = cat(1,post{:,c});
        tdesign{c} = cat(1,design{:,c});
      end
    
      p = cell(1,length(tpost));
      for c=1:length(tpost)
    
        if ~isfield(options,'test') || strcmp(options.test,'binomial_test')
        % one-sided binomial test with automatic bonferroni correction
        
          if obj.verbose
            fprintf('performing one-sided biniomial test with bonferroni correction\n');
          end
          
          % compute class with highest prior probability
          nclasses = size(tpost{c},2);
          priors = zeros(1,nclasses);
          for k=1:nclasses
            priors(k) = priors(k) + sum(tdesign{c}(:,1)==k);
          end
          [mxx,clss] = max(priors);
          
          rndpost = zeros(size(tpost{c}));
          rndpost(:,clss) = 1;
          
          [r,p{c},level] = validator.binomial_test(tpost{c},tdesign{c},rndpost,tdesign{c},'twosided',false,'bonferroni',length(tpost));
          
          if r
            fprintf('dataset %d: null hypothesis rejected (%g<%g);\nsignificant difference from ',c,p{c},level);
          else
            fprintf('dataset %d: null hypothesis not rejected (%g>%g);\nno significant difference from ',c,p{c},level);
          end
          
          if (nargin == 2 && random)
            fprintf('random classification.\n');
          else
            fprintf('majority classification (class %d with prior of %g)\n',clss,mxx/sum(priors));
          end
          
        end
      end
      
      if length(p) == 1
        p = p{1};
      end
      
    end
    
  end
  
  methods(Access = protected)
 
    
    function p = getpredictor(obj)
       % get the predictor (ending) mvmethod

       if iscell(obj.procedure)
         if iscell(obj.procedure{1}.mvmethods{end})
           p = obj.procedure{1}.mvmethods{end}{1};
         else
           p = obj.procedure{1}.mvmethods{end};
         end
       else
         if iscell(obj.procedure.mvmethods{end})
           p = obj.procedure.mvmethods{end}{1};
         else
           p = obj.procedure.mvmethods{end};
         end
       end
      
    end
    
  end
  
  methods(Static)
    
    
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
