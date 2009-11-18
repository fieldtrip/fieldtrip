classdef l1lr < classifier
%L1 regularized logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: l1lr.m,v $
%

    properties

        model; 
        
        lambda = 1;
        
    end

    methods
       function obj = l1lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
         
         [data,design] = obj.check_input(data,design);
         
         if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
         if obj.nclasses ~= 2, error('L2LR expects binary class labels'); end

         data = [ones(size(data,1),1) data];

         design = design - 1; % convert to -1 / + 1
         
         funObj = @(w)LogisticLoss(w,data,design);
         w_init = zeros(size(data,2),1);
         
         if obj.verbose
           fprintf('\nComputing L1-Regularized Logistic Regression Coefficients...\n');
         end
           
         lambdas = obj.lambda*ones(size(data,2),1);
         lambdas(1) = 0; % Do not penalize bias variable
         
         obj.model = L1GeneralProjection(funObj,w_init,lambdas);


       end

       function post = test(obj,data)

         data = obj.check_input(data);
         
         data = [ones(size(data,1),1) data];
         
         post = 1 ./ (1 + exp(data * obj.model));
         post = [post 1 - post];
         
       end              
       
       function m = getmodel(obj,label,dims)
         % return the parameters wrt a class label in some shape 
         
           m = obj.model(2:end)'; % ignore bias term
         
         if nargin == 3 && numel(m) == prod(dims)
           m = reshape(m,dims);
         end
         
       end
       
    end
end 
