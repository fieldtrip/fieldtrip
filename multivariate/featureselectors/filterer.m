classdef filterer < featureselector
%FILTERER filtering approach to feature selection
%
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
%   PARAMETERS:
%    criterion; % keeps track of the evaluation metric
%    value; % filter function outputs
%    order; % ordering of the variables
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
  
    validator = []; % e.g., crossvalidator('procedure',mva({nb()}),'cvfolds',0.9);
    metric = 'accuracy'; % evaluation metric
   
    nfeatures = Inf; % maximum number of used features
        
    % filter is used to order the features
    %
    % options:
    % [10 100 2 ...] : uses the ordering specified
    % @meandiff : order according to the norm of differences of means
    % @mutual_information : order according to the mutual information (default)
    % @anova : order according to anova test    
    filter = @mutual_information;
    
    
  end
  
  methods
    function obj = filterer(varargin)
      
      obj = obj@featureselector(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
      
      maxfeatures = min(obj.nfeatures,size(X,2));
      
      sfilt = func2str(obj.filter);
      
      if obj.verbose
        if isempty(obj.validator)
          fprintf('selecting %d features based on %s filter\n',maxfeatures,sfilt);
        else
          fprintf('selecting optimal features based on %s filter\n',sfilt);
        end
      end
        
      % compute function values
      nfeatures = size(X,2);
      p.value = zeros(1,nfeatures);
      for j=1:nfeatures
        
        fprintf('computing %s filter for feature %d of %d\n',sfilt,j,nfeatures);
        
        p.value(j) = obj.filter(X(:,j),Y);
      end
      
      % get ordering of the features
      [tmp,p.order] = sort(p.value,'descend');
      
      if isempty(obj.validator)
        % select the best n features
        p.subset = p.order(1:maxfeatures);
      else
        % use procedure to determine the best feature subset
        [p.subset,p.criterion] = obj.select_features(X,Y,p.order);
      end
      
    end
    
    function [subset,criterion] = select_features(obj,X,Y,features)
      % select features
      % assumes availability of evaluation function and feature
      % ordering
      
      subset = [];
      metric = -Inf;
            
      criterion = zeros(1,length(features));
      for f=1:min(obj.nfeatures,length(features))
        
        if obj.verbose
          fprintf('evaluating %d out of %d features; ',f,length(features));
        end
        
        cv = obj.validator.validate(X(:,features(1:f)),Y);
        
        m = evaluate(cv.Y,cv.Y,'metric',obj.metric);
        
        if obj.verbose, fprintf('criterion: %.2g\n',m); end
        
        if m > metric
          subset = features(1:f);
          metric = m;
        end
        
        criterion(f) = m;
        
      end
      
    end
    
  end
  
  
end

