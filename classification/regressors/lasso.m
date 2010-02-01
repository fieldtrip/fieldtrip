classdef lasso < regressor
%LASSO using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven

    properties
        
      lambda = 100;
      
    end

    methods
      
       function obj = lasso(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       
       function p = estimate(obj,data,design)
     
         data = [data.X ones(data.nsamples,1)];
         
         lambdas = obj.lambda*ones(size(data,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(data'*data + diag(lambdas));

         funObj = @(w)GaussianLoss(w,data,design.X); % Loss function that L1 regularization is applied to
         
         w_init = R\(R'\(data'*design.X)); % Initial value for iterative optimizer
         
         if obj.verbose
           fprintf('\nComputing LASSO Coefficients...\n');
         end
         
         p.model = L1GeneralProjection(funObj,w_init,lambdas);

       end
       
       function post = map(obj,data)       
       
         post = dataset([data.X ones(data.nsamples,1)] * obj.params.model);

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters 
                    
         m = {full(obj.params.model(:,1:(end-1)))}; % ignore bias term
         desc = {'unknown'};
         
       end

    end
end 
