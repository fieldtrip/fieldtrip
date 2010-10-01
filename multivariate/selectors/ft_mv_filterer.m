classdef ft_mv_filterer < ft_mv_selector
%FT_MV_FILTERER filtering approach to feature selection
%
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
%   EXAMPLE:
%   ft_mv_filterer('validator',ft_mv_crossvalidator('procedure',ft_mv_naive,'metric','accuracy'))
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    maxfeatures = Inf; % maximum number of used features
        
    filter = @(x,y)(abs(corr(x,y))); % the used filter
    
    value % the value  of the filter for each feature
    
    order % the ordering of the features according to value
    
    criterion % the value of the validator for feature 1:N    
    
  end
  
  methods
    
    function obj = ft_mv_filterer(varargin)
      
      obj = obj@ft_mv_selector(varargin{:});
      
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
        for j=1:nf
          
          if obj.verbose
            fprintf('computing %s filter for feature %d of %d\n',obj.filter,j,nf);
          end
          
          obj.value(j) = obj.filter(X(:,j),Y);
          
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
        for j=1:nf
          
          if obj.verbose
            fprintf('computing %s filter for feature %d of %d\n',obj.filter,j,nf);
          end
          
          obj.value(j) = obj.filter(X(:,j),Y);
          
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
          
          cv = obj.validator.validate(X(:,obj.order(1:f)),Y);
          
          m = cv.evaluate;
          
          if obj.verbose, fprintf('criterion: %.2g\n',m); end
          
          if m > metric
            obj.subset = obj.order(1:f);
            metric = m;
          end
          
          obj.criterion(f) = m;
          
        end
      end
      
    end
   
    
    function [model,desc] = model(obj)
      % return logical array with ones indicating selected features
      
      indim = obj.indims;
      if length(indim) == 1, indim = [1 indim]; end
      
      m = zeros(indim);
      m(obj.params.subset) = 1;      
      model{1} = m;      
      desc{1} = {'logical array with ones indicating selected features'};
      
      m = zeros(indim);
      m(:) = obj.params.value;
      model{2} = m;
      desc{2} = {sprintf('importance according to filter %s',obj.filter)};
      
    end
    
  end
  
end