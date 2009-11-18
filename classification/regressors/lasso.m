classdef lasso < regressor
%LASSO using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: lasso.m,v $
%

    properties
        
      model     
      
      lambda = 100;
      
    end

    methods
       function obj = lasso(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function obj = train(obj,data,design)
        
         % add bias term
         data = [data ones(size(data,1),1)];
         
         lambdas = obj.lambda*ones(size(data,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(data'*data + diag(lambdas));

         funObj = @(w)GaussianLoss(w,data,design); % Loss function that L1 regularization is applied to
         
         w_init = R\(R'\(data'*design)); % Initial value for iterative optimizer
         
         if obj.verbose
           fprintf('\nComputing LASSO Coefficients...\n');
         end
         
         obj.model = L1GeneralProjection(funObj,w_init,lambdas);

       end
       
       function post = test(obj,data)       
       
         % add bias term
         data = [data ones(size(data,1),1)];
       
         post = data * obj.model;

       end
       
       function m = getmodel(obj,label,dims)
         % return the parameters wrt a class label in some shape 
                    
         m = full(obj.model(:,1:(end-1))); % ignore bias term

         if nargin == 3 && numel(m) == prod(dims)
           m = reshape(m,dims);
         end
         
       end

    end
end 
