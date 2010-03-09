classdef elasticlr < classifier
%ELASTICLR elastic net logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

      lambdaL1 = 1;
      lambdaL2 = 1;
      
    end

    methods
      
       function obj = elasticlr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = estimate(obj,X,Y)
               
         if obj.nunique(Y) ~= 2, error('elasticLR expects binary class labels'); end

         X = [ones(size(X,1),1) X];

         design = 2*Y - 3; % convert to -1 / + 1
         
         funObj = @(w)LogisticLoss(w,X,design);
         
         if ~isempty(obj.params) && isfield(obj.params,'model') && ~isempty(obj.params.model) 
             % initial model has been specified
             w_init = obj.params.model;
         else
             w_init = zeros(size(X,2),1);
         end
         
         if obj.verbose
           fprintf('\nComputing Elastic-Net Logistic Regression Coefficients...\n');
         end
           
         lambdasL1 = obj.lambdaL1*ones(size(X,2),1);
         lambdasL1(1) = 0; % Do not penalize bias variable
         
         lambdasL2 = obj.lambdaL2*ones(size(X,2),1);
         lambdasL2(1) = 0; % Do not penalize bias variable

         funObjL2 = @(w)penalizedL2(w,funObj,lambdasL2);

         opt.verbose = 0;
         p.model = L1GeneralProjection(funObjL2,w_init,lambdasL1,opt);

       end

       function Y = map(obj,X)

         Y = 1 ./ (1 + exp([ones(size(X,1),1) X] * obj.params.model));
         Y = [Y 1 - Y];  
          
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
         m = {obj.params.model(2:end)'}; % ignore bias term
         desc = {'unknown'};
         
       end
       
    end
end 
