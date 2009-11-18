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
        
        [res all] = evaluate(obj.post,obj.design,varargin{:});
        
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

        [result all] = evaluate(obj.post,obj.design,varargin{:});
        
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
      [r,p,level] = significance(tpost,tdesign,rndpost,tdesign,'twosided',false);
      
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
    
    
  end 
    
end
