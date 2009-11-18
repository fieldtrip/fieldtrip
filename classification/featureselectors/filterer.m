classdef filterer < featureselector
%FILTERER filtering approach to feature selection
%
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
% OPTIONS
%
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: filterer.m,v $
%

    properties
        
        % the filter function should accept a dataset and design and output
        % a measure of desirability for that dataset. In principle, a
        % classifier could be used to compute univariate accuracies as a
        % measure (also see wrapper). 
        %
        % We can use any of:
        %
        % [] : uses the validator to order the features (not yet implemented)
        % [10 100 2 ...] : uses the ordering specified
        % @meandiff : order according to the norm of differences of means
        % @mutual_information : order according to the mutual information
        % @anova : order according to anova test
        
        filter = @meandiff
    
        value; % filter function outputs
        order; % ordering of the variables
        top; % selection of the best features (false by default)
     end

    methods
        function obj = filterer(varargin)

            obj = obj@featureselector(varargin{:});
                    
        end
        function obj = train(obj,data,design)
            
            if obj.verbose
                if obj.top
                    fprintf('selecting %d features based on %s filter\n',obj.top,func2str(obj.filter));
                else
                    fprintf('selecting optimal features based on %s filter\n',func2str(obj.filter));
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
                    obj.value(j) = obj.filter(data(:,j),design);
                end
                
                % get ordering of the features
                [a,obj.order] = sort(obj.value,'descend');

                if obj.top
                    % select the best top features
                    obj.subset = obj.order(1:obj.top);
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
                for f=1:min(obj.maxfeatures,length(features))
                   
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

