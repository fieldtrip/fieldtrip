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
               
         if design.nunique ~= 2, error('elasticLR expects binary class labels'); end

         design = design.X;
         
         X = [ones(data.nsamples,1) data.X];

         design = 2*design - 3; % convert to -1 / + 1
         
         funObj = @(w)LogisticLoss(w,X,design);
         w_init = zeros(data.nfeatures+1,1);

         if obj.verbose
           fprintf('\nComputing Elastic-Net Logistic Regression Coefficients...\n');
         end
           
         lambdasL1 = obj.lambdaL1*ones(data.nfeatures+1,1);
         lambdasL1(1) = 0; % Do not penalize bias variable
         
         lambdasL2 = obj.lambdaL2*ones(data.nfeatures+1,1);
         lambdasL2(1) = 0; % Do not penalize bias variable

         funObjL2 = @(w)penalizedL2(w,funObj,lambdasL2);

         obj.model = L1GeneralProjection(funObjL2,w_init,lambdasL1,'verbose',0);

       end

       function post = test(obj,data)

         post = 1 ./ (1 + exp([ones(data.nsamples,1) data.X] * obj.model));
         post = dataset([post 1 - post]);  
          
       end              
       
       function m = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
         m = obj.model(1:(end-1))'; % ignore bias term
         
       end
       
    end
end 
