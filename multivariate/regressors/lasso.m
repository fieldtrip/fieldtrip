classdef lasso < regressor
%LASSO using Mark Schmidt's L1General package
%
% NOTE: 
% using elasticnet with lambda = 0 and nu ~= 0 can be faster
%
%   Copyright (c) 2009, Marcel van Gerven

    properties
        
      lambda = 100;
      
    end

    methods
      
       function obj = lasso(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       
       function p = estimate(obj,X,Y)
     
         X = [X ones(size(X,1),1)];
         
         lambdas = obj.lambda*ones(size(X,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(X'*X + diag(lambdas));

         funObj = @(w)GaussianLoss(w,X,Y); % Loss function that L1 regularization is applied to
         
         w_init = R\(R'\(X'*Y)); % Initial value for iterative optimizer
         
         if obj.verbose
           fprintf('\nComputing LASSO Coefficients...\n');
         end
         
         p.model = L1GeneralProjection(funObj,w_init,lambdas);

       end
       
       function Y = map(obj,X)       
       
         Y = [X ones(size(X,1),1)] * obj.params.model;

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters 
                    
         m = {full(obj.params.model(:,1:(end-1)))}; % ignore bias term
         desc = {'unknown'};
         
       end

    end
end 
