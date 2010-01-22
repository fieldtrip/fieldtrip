classdef filterer < featureselector
%FILTERER filtering approach to feature selection
%
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: filterer.m,v $
%

  properties
  
    validator = []; % e.g., crossvalidator('procedure',clfproc({nb()}),'cvfolds',0.9);
    metric = 'accuracy'; % evaluation metric
    criterion; % keeps track of the evaluation metric
    
    nfeatures = Inf; % maximum number of used features
        
    % filter is used to order the features
    %
    % options:
    % [10 100 2 ...] : uses the ordering specified
    % @meandiff : order according to the norm of differences of means
    % @mutual_information : order according to the mutual information (default)
    % @anova : order according to anova test    
    filter = @mutual_information;
    
    value; % filter function outputs
    order; % ordering of the variables
    
  end
  
  methods
    function obj = filterer(varargin)
      
      obj = obj@featureselector(varargin{:});
      
    end
    function obj = train(obj,data,design)
      
      obj.nfeatures = min(obj.nfeatures,data.nfeatures);
      
      data = data.collapse();
      design = design.collapse();
      
      sfilt = func2str(obj.filter);
      
      if obj.verbose
        if isempty(obj.validator)
          fprintf('selecting %d features based on %s filter\n',obj.nfeatures,sfilt);
        else
          fprintf('selecting optimal features based on %s filter\n',sfilt);
        end
      end
      
      if iscell(data)
        
        cvalue = cell(1,length(data));
        corder = cell(1,length(data));
        ccriterion = cell(1,length(data));
        csubset = cell(1,length(data));
        
        for c=1:length(data)
          obj = obj.train(data{c},design{c});
          cvalue{c} = obj.value;
          corder{c} = obj.order;
          ccriterion{c} = obj.criterion;
          csubset{c} = obj.subset;
        end
        
        obj.value = cvalue;
        obj.order = corder;
        obj.criterion = ccriterion;
        obj.subset = csubset;
        
      else
        
        % compute function values
        nfeatures = size(data,2);
        obj.value = zeros(1,nfeatures);
        for j=1:nfeatures
          
          fprintf('computing %s filter for feature %d of %d\n',sfilt,j,nfeatures);
          
          obj.value(j) = obj.filter(data(:,j),design);
        end
        
        % get ordering of the features
        [a,obj.order] = sort(obj.value,'descend');
        
        if isempty(obj.validator)
          % select the best n features
          obj.subset = obj.order(1:obj.nfeatures);
        else
          % use procedure to determine the best feature subset
          [obj.subset,obj.criterion] = obj.select_features(data,design);
        end
      end
      
    end
    
    function [subset,criterion] = select_features(obj,data,design)
      % select features
      % assumes availability of evaluation function and feature
      % ordering
      
      subset = [];
      metric = -Inf;
      
      % use a more fancy procedure
      features = obj.order;
      
      criterion = zeros(1,length(features));
      for f=1:min(obj.nfeatures,length(features))
        
        if obj.verbose
          fprintf('evaluating %d out of %d features; ',f,length(features));
        end
        
        cv = obj.validator.validate(data(:,features(1:f)),design);
        
        m = evaluate(cv.post,cv.design,'metric',obj.metric);
        
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

