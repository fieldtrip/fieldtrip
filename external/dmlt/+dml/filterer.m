classdef filterer < dml.method
% FILTERER filtering approach to feature selection.
%
%   DESCRIPTION
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
%   If filter is a string, e.g. 'ttest' then it calls rankfeatures in the bioinformatics
%   toolbox. If it is a function it will just evaluate that function on
%   each feature and the output.
%
%   EXAMPLE
%   dml.filterer('validator',dml.crossvalidator('mva',dml.naive,'stat','accuracy'))
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    maxfeatures = Inf; % maximum number of used features
        
    filter = @(x,y)(abs(corr(x,y))); % the used filter
    
    value % the value  of the filter for each feature
    
    order % the ordering of the features according to value
    
    criterion % the value of the validator for feature 1:N    
    
    validator % the used crossvalidator
    
    subset % subset of selected features
    
  end
  
  methods
    
    function obj = filterer(varargin)
      
      obj = obj@dml.method(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
            
      maxf = min(obj.maxfeatures,size(X,2));
      
      if obj.verbose
        if isempty(obj.validator)
          fprintf('selecting %d features based on %s filter\n',maxf,obj.filter);
        else
          fprintf('selecting optimal features based on %s filter\n',obj.filter);
        end
      end      
      
      if isempty(obj.validator)
      % just take the best N features based on a univariate criterion
      
        % compute function values
        nf = size(X,2);
        obj.value = zeros(1,nf);
        if ischar(obj.filter)
          
          if obj.verbose
            fprintf('computing %s filter using rankfeatures\n');
          end
          Z = rankfeatures(X',Y','Criterion',obj.filter);
          for j=1:nf
            obj.value(j) = Z(j);
          end
          
        else
          
          for j=1:nf
            
            if obj.verbose
              fprintf('computing %s filter for feature %d of %d\n',obj.filter,j,nf);
            end
            
            obj.value(j) = obj.filter(X(:,j),Y);
          end
          
        end
        
        % get ordering of the features
        [tmp,obj.order] = sort(obj.value,'descend');
      
        % select the best n features
        obj.subset = obj.order(1:maxf);
      
      else
        % currently we use a simple approach. A better approach would be
        % the following:
        %
        % use procedure to determine the best feature subset.
        % first we use the procedure to compute the univariate filter
        % values per fold. Then we determine the average best number of features
        % over folds using the crossvalidation procedure. Then we recompute the 
        % univariate filter values for all data and return the optimal
        % subset
        
        % compute function values
        nf = size(X,2);
        obj.value = zeros(1,nf);
        if ischar(obj.filter)
          
          if obj.verbose
            fprintf('computing %s filter using rankfeatures\n');
          end
          Z = rankfeatures(X',Y','Criterion',obj.filter);
          for j=1:nf
            obj.value(j) = Z(j);
          end
          
        else
          for j=1:nf
            
            if obj.verbose
              fprintf('computing %s filter for feature %d of %d\n',obj.filter,j,nf);
            end
            
            obj.value(j) = obj.filter(X(:,j),Y);
          
          end
          
        end
        
        % get ordering of the features
        [tmp,obj.order] = sort(obj.value,'descend');
        
        obj.subset = [];
        metric = -Inf;
        
        obj.criterion = zeros(1,maxf);
        for f=1:maxf
          
          if obj.verbose
            fprintf('evaluating %d out of %d features; ',f,maxf);
          end
          
          cv = obj.validator.train(X(:,obj.order(1:f)),Y);
          
          m = cv.statistic;
          
          if obj.verbose, fprintf('criterion: %.2g\n',m); end
          
          if m > metric
            obj.subset = obj.order(1:f);
            metric = m;
          end
          
          obj.criterion(f) = m;
          
        end
      end
      
    end
   
    function Y = test(obj,X)
            
      Y = X(:,obj.subset);
      
    end
    
    function m = model(obj)
      % return logical array with ones indicating selected features
      % m.selected logical array with ones indicating selected features
      % m.importance importance according to the used filter
      
      indim = obj.indims;
      if length(indim) == 1, indim = [1 indim]; end
      
      m1 = zeros(indim);
      m1(obj.subset) = 1;      
      m.selected = m1;    
      
      m1 = zeros(indim);
      m1(:) = obj.value;
      m.importance = m1;
      
    end
    
  end
  
end

