classdef elasticnet < regressor
%ELASTICNET using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: elasticnet.m,v $
%

    properties
        
      model     
      
      lambdaL1 = 100;
      lambdaL2 = 100;

    end

    methods
       function obj = elasticnet(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function obj = train(obj,data,design)
        
         data = data.collapse();
         design = design.collapse();
         
         % add bias term
         data = [data ones(size(data,1),1)];
        
         lambdasL2 = obj.lambdaL2*ones(size(data,2),1);
         lambdasL1 = obj.lambdaL1*ones(size(data,2),1);

         funObj = @(w)GaussianLoss(w,data,design); % Loss function that L1 regularization is applied to
         
         R = chol(data'*data + diag(lambdasL1));
         w_init = R\(R'\(data'*design)); % Initial value for iterative optimizer
         
         penalizedFunObj = @(w)penalizedL2(w,funObj,lambdasL2);

         if obj.verbose
           fprintf('\nComputing Elastic Net Coefficients...\n');
         end
         
         obj.model = L1GeneralProjection(penalizedFunObj,w_init,lambdasL1);

       end
       
       function post = test(obj,data)       
       
         % add bias term
         data = [data.collapse() ones(data.nsamples,1)];       
         post = dataset(data * obj.model);

       end
       
       function m = getmodel(obj)
         % return the parameters wrt a class label in some shape 
                    
         m = full(obj.model(1:(end-1))); % ignore bias term
         
       end

    end
end 
