classdef elasticlr < classifier
%ELASTICLR elastic net logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: elasticlr.m,v $
%

    properties

        model; 
        
        lambdaL1 = 1;
        lambdaL2 = 1;
      
    end

    methods
       function obj = elasticlr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
         
         [data,design] = obj.check_data(data,design);
         
         if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
         if obj.nclasses ~= 2, error('L2LR expects binary class labels'); end

         data = [ones(size(data,1),1) data];

         design = design - 1; % convert to -1 / + 1
         
         funObj = @(w)LogisticLoss(w,data,design);
         w_init = zeros(size(data,2),1);

         if obj.verbose
           fprintf('\nComputing Elastic-Net Logistic Regression Coefficients...\n');
         end
           
         lambdasL1 = obj.lambdaL1*ones(size(data,2),1);
         lambdasL1(1) = 0; % Do not penalize bias variable
         lambdasL2 = obj.lambdaL2*ones(size(data,2),1);
         lambdasL2(1) = 0; % Do not penalize bias variable

         funObjL2 = @(w)penalizedL2(w,funObj,lambdasL2);

         obj.model = L1GeneralProjection(funObjL2,w_init,lambdasL1);

       end

       function post = test(obj,data)

         data = obj.check_input(data);
         
         data = [ones(size(data,1),1) data];
         
         post = 1 ./ (1 + exp(data * obj.model));
         post = [post 1 - post];
          
       end              
       
       function m = getmodel(obj,label,dims)
         % return the parameters wrt a class label in some shape 
         
           m = obj.model(1:(end-1))'; % ignore bias term
         
         if nargin == 3 && numel(m) == prod(dims)
           m = reshape(m,dims);
         end
         
       end
       
    end
end 
