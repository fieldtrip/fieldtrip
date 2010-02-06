classdef regressor < predictor
%REGRESSOR abstract regression method class
%
% A regressor takes a variable number of arguments upon construction. 
% During operation, the regressor takes data and
% produces predicted responses as an N x 1 matrix for N examples.
% 
% Subclasses should implement the train and test functions.
%
% OPTIONS
%   'nclasses'  : always nan for a regressor
%   'nfeatures' : number of features (determined from data)
%   'nexamples' : number of examples (determined from data)
%
% SEE ALSO
%   doc regressors
%
%   Copyright (c) 2008, Marcel van Gerven
     
    
    methods
        
        function obj = regressor(varargin)
               
          % parse options
          for i=1:2:length(varargin)
            if ismember(varargin{i},fieldnames(obj))
              obj.(varargin{i}) = varargin{i+1};
            end
          end

        end
        
        function reg = predict(obj,data)
           % convert posterior dataset into predictions (mean value)
           
           reg = obj.test(data).X;
           reg = reg(:,1);
           
        end
    end
    
    methods(Static)
      
      function metric = evaluate(post,design,varargin)
        %EVALUATE evaluation metrics for regressors
        %
        %   metric = evaluate(post,design,varargin)
        %
        %   parameter 'metric' determines the evaluation criterion:
        %
        %       'rss' : residual sum of squares
        %       'mse' : mean squared error
        %       'circrss' : circular residual sum of squares
      
         options = varargin2struct(varargin);
      
         if ~isfield(options,'metric'), options.metric = 'accuracy'; end
         
         if isa(post,'dataset')
           metric = compute_metric(post.X(:,1),design.X(:,1),options);
         else
           metric = compute_metric(post(:,1),design(:,1),options);
         end
      
         function met = compute_metric(post,design,cfg)
      
           met = [];
           switch lower(cfg.metric)
             
             case 'rss' % residual sum of squares
               
               met = sum((design - post).^2);
             
             case 'mse' % mean squared error
               
               met = sum((design - post).^2);
              
             case 'circrss' % circular residual sum of squares; assumes data is in the range -pi .. pi
               
               met = mod(mod(design,2*pi) - mod(post,2*pi),2*pi);
               met = min(met,2*pi - met);
               met = sum(met.^2);
           
             otherwise
                error(['unsupported option ',cfg.metric]);
           
           end
           
         end
      end
    end
      
end
