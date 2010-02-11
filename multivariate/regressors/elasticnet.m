classdef elasticnet < regressor
%ELASTICNET using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven


    properties
        
      lambdaL1 = 100;
      lambdaL2 = 100;

    end

    methods
      
       function obj = elasticnet(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       
       function p = estimate(obj,X,Y)
        
         X = [X ones(size(X,1),1)];
        
         lambdasL2 = obj.lambdaL2*ones(size(X,2),1);
         lambdasL1 = obj.lambdaL1*ones(size(X,2),1);

         funObj = @(w)GaussianLoss(w,X,Y); % Loss function that L1 regularization is applied to
         
         R = chol(X'*X + diag(lambdasL1));
         w_init = R\(R'\(X'*Y)); % Initial value for iterative optimizer
         
         penalizedFunObj = @(w)penalizedL2(w,funObj,lambdasL2);

         if obj.verbose
           fprintf('\nComputing Elastic Net Coefficients...\n');
         end
         
         p.model = L1GeneralProjection(penalizedFunObj,w_init,lambdasL1);

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
