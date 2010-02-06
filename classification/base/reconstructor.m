classdef reconstructor < predictor
%RECONSTRUCTOR reconstructor class
%
% a reconstructor takes multidimensional inputs and outputs and
% reconstructs the outputs. This can be done by simply applying classifiers
% to individual outputs but also by much more sophisticated methods
% 
% Subclasses should implement the train and test functions and optionally
% the getmodel function which reshapes parameters to something
% interpretable
%
% Copyright (c) 2010, Marcel van Gerven
  

  methods
        
        function obj = reconstructor(varargin)
     
          % parse options
          for i=1:2:length(varargin)
            if ismember(varargin{i},fieldnames(obj))
              obj.(varargin{i}) = varargin{i+1};
            end
          end

        end        
        
        function rec = predict(obj,data)
          % convert outcome dataset into predictions
          
          rec = obj.test(data).X; 
        end
        
  end
    
  methods(Static)
      
      function metric = evaluate(post,design,varargin)
        %EVALUATE evaluation metrics for reconstructors
        %
        %   metric = evaluate(post,design,varargin)
        %
        %   parameter 'metric' determines the evaluation criterion:
        %
        %       'rss' : residual sum of squares
        %       'mse' : mean squared error
        %       'correlation' : correlation between real and predicted reconstructions
      
         options = varargin2struct(varargin);
      
         if ~isfield(options,'metric'), options.metric = 'accuracy'; end
         
         if isa(post,'dataset')
           metric = compute_metric(post.X,design.X,options);
         else
           metric = compute_metric(post,design,options);
         end
      
         function met = compute_metric(post,design,cfg)
      
           met = [];
           switch lower(cfg.metric)
             
             case 'rss' % residual sum of squares
               
               met = sum((design - post).^2);
             
             case 'mse' % mean squared error
               
               met = sum((design - post).^2);
              
             case 'correlation' % correlation
               
               met = corrcoef(post(:),design(:));
               met = met(2);
               
             case 'cityblock'
               
             case 'cosine'
               
               
             otherwise
                error(['unsupported option ',cfg.metric]);
           
           end
           
         end
      end
    end
    
end
