classdef elasticnet < regressor
%ELASTICNET regressor
%
%   Copyright (c) 2010, Tom Heskes, Marcel van Gerven


    properties
        
      lambda = 0.4; % ridge penalty
      nu = 0.01; % L1 penalty

    end

    methods
      
       function obj = elasticnet(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       
       function p = estimate(obj,X,Y)        
         
         [beta,beta0] = elastic(X',Y',obj.nu,obj.lambda);
         
         p.model = [beta; beta0];

         if ~any(p.model)
           warning('elasticnet returned zero vector as model; too strong regularization?');
         end
         
       end
       
       function Y = map(obj,X)       
       
         Y = [X ones(size(X,1),1)] * obj.params.model;

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
                    
         m = {obj.params.model(1:(end-1))'}; % ignore bias term
         desc = {'unknown'};
                    
       end

    end
end 
