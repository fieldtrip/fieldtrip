classdef l2lr < classifier
%L2 regularized logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: l2lr.m,v $
%

    properties

        model; 
        
        lambda = 1;
        
    end

    methods
       function obj = l2lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)

        if design.nunique > 2, error('L2LR expects binary class labels'); end

         if obj.verbose
           fprintf('\nComputing L2-Regularized Logistic Regression Coefficients...\n');
         end
         
         design = design.collapse();
               
         % make foolproof
         if all(design(:,1) == 1)
           obj.model = -inf;
         elseif all(design(:,1) == 2)
           obj.model = inf;
         else
           
           X = [ones(data.nsamples,1) data.collapse()];
           
           design = 2*design - 3; % convert to -1 / + 1
           
           funObj = @(w)LogisticLoss(w,X,design);
           w_init = zeros(data.nfeatures+1,1);
           
           lambdas = obj.lambda*ones(data.nfeatures+1,1);
           lambdas(1) = 0; % Do not penalize bias variable
           
           funObjL2 = @(w)penalizedL2(w,funObj,lambdas);
           
           mfOptions.Method = 'lbfgs';
           mfOptions.display = 0;
           obj.model = minFunc(funObjL2,w_init,mfOptions);
       
         end
         
       end

       function post = test(obj,data)

        if isinf(obj.model(1)) 
           if obj.model < 0
            post = dataset([ones(data.nsamples,1) zeros(data.nsamples,1)]);
           else
             post = dataset([zeros(data.nsamples,1) ones(data.nsamples,1)]);
           end
         else
           post = 1 ./ (1 + exp([ones(data.nsamples,1) data.collapse()] * obj.model));
           post = dataset([post 1 - post]);
        end
         
       end              
       
       function m = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = obj.model(1:(end-1))'; % ignore bias term
         
       end
       
    end
end 
