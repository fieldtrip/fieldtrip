classdef l1lr < classifier
%L1 regularized logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

        lambda = 1;
        
    end

    methods
      
       function obj = l1lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = estimate(obj,X,Y)
                 
         if obj.nunique(Y) ~= 2, error('L1LR expects binary class labels'); end

         if obj.verbose
           fprintf('\nComputing L1-Regularized Logistic Regression Coefficients...\n');
         end
         
         % foolproof
         if all(Y(:,1) == 1)
           p.model = -inf;
         elseif all(Y(:,1) == 2)
           p.model = inf;
         else
           
           X = [ones(size(X,1),1) X];
           
           Y = 2*Y - 3; % convert to -1 / + 1
           
           lambdas = obj.lambda*ones(size(X,2),1);
           lambdas(1) = 0; % Do not penalize bias variable
           
           funObj = @(w)LogisticLoss(w,X,Y);
           w_init = zeros(size(X,2),1);
            
           opt.verbose = 0;
           p.model = L1GeneralProjection(funObj,w_init,lambdas,opt);
         end

       end

       function Y = map(obj,X)

         if isinf(obj.params.model(1)) 
           if obj.params.model < 0
             Y = [ones(size(X,1),1) zeros(size(X,1),1)];
           else
             Y = [zeros(size(X,1),1) ones(size(X,1),1)];
           end
         else
           Y = 1 ./ (1 + exp([ones(size(X,1),1) X] * obj.params.model));
           Y = [Y 1 - Y];
         end
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.params.model(2:end)'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
