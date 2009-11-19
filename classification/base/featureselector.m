classdef featureselector < clfmethod
%FEATURESELECTOR abstract featureselector method class
%
% A featureselector is a handle class that takes a variable number of arguments 
% upon construction. During operation, the featureselector takes data and
% produces a reduced dataset with (typically) a smaller number of features.
% A featureselector differs from a preprocessor in the sense that it
% typically requires a VALIDATOR and an evaluation CRITERION in order to
% determine which features to use.
%
% OPTIONS
%   'subset'    : the subset which can be set manually using this class
%   'validator' : validation procedure to use
%   'verbose'   : output comment if true
%   'metric'    : evaluation metric
%
% Subclasses should implement the train and test functions and possibly
% getmodel.
%
% SEE ALSO
% doc featureselectors
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: featureselector.m,v $
%
    properties
        subset = []; % the feature subset that is to be selected
        validator = []; % e.g., crossvalidator('procedure',clfproc({nb()}),'cvfolds',0.9);
        metric = 'accuracy'; % evaluation metric
        criterion; % keeps track of the evaluation metric
        maxfeatures = Inf; % maximum number of used features   
    end

    methods
        
        function obj = featureselector(varargin)

            % parse options 
            for i=1:2:length(varargin)
              if ismember(varargin{i},fieldnames(obj))
                obj.(varargin{i}) = varargin{i+1};
              end
            end
            
            assert(isempty(obj.validator) || isa(obj.validator,'validator'));
        
        end

        function obj = train(obj,data,design)
            
        end
        
        function data = test(obj,data)
            
            if iscell(data)
                
                for c=1:length(data)

                    % return data for the subset
                    data{c} = data{c}(:,obj.subset{c});
                end
            else
                data = data(:,obj.subset);
            end
        end
        
    end
    
end
