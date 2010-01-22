classdef featureselector < clfmethod
%FEATURESELECTOR featureselector method class
%
% During operation, the featureselector takes data and
% produces a reduced dataset with (typically) a smaller number of features.
%
% OPTIONS
%   'subset'    : the subset which can be set manually using this class
%
% Subclasses should implement the train and/or test functions and possibly
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
      
      subset = []; % the feature subset that is to be used
      
    end

    methods
      
      function obj = featureselector(varargin)
        
        % parse options
        for i=1:2:length(varargin)
          if ismember(varargin{i},fieldnames(obj))
            obj.(varargin{i}) = varargin{i+1};
          end
        end
              
      end
      
      function obj = train(obj,data,design)
        
      end
      
      function data = test(obj,data)
        
        if obj.verbose
          fprintf('using %d out of %d features\n',numel(obj.subset),data.nfeatures);
        end
        
        if iscell(data)
          
          for c=1:length(data)
            
            % return data for the subset
            data{c} = data{c}.collapse();
            data{c} = dataset(data{c}(:,obj.subset{c}));
          end
        else
          data = data.collapse();
          data = dataset(data(:,obj.subset));
        end
        
      end
      
      function model = getmodel(obj,label)
        % return used subset
        
        model = obj.subset;
        
      end
      
    end
    
end
