classdef l2lr < classifier
%L2 regularized logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

        lambda = 1;
        
    end

    methods
       function obj = l2lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function p = estimate(obj,X,Y)

        if obj.nunique(Y) > 2, error('L2LR expects binary class labels'); end

         if obj.verbose
           fprintf('\nComputing L2-Regularized Logistic Regression Coefficients...\n');
         end
               
         % make foolproof
         if all(Y(:,1) == 1)
           p.model = -inf;
         elseif all(Y(:,1) == 2)
           p.model = inf;
         else
           
           X = [ones(size(X,1),1) X];
           
           Y = 2*Y - 3; % convert to -1 / + 1
           
           funObj = @(w)LogisticLoss(w,X,Y);
           
           if ~isempty(obj.params) && isfield(obj.params,'model') && ~isempty(obj.params.model) 
             % initial model has been specified
             w_init = obj.params.model;
           else
             w_init = zeros(size(X,2),1);
           end
           
           lambdas = obj.lambda*ones(size(X,2),1);
           lambdas(1) = 0; % Do not penalize bias variable
           
           funObjL2 = @(w)penalizedL2(w,funObj,lambdas);
           
           opt.Method = 'lbfgs';
           opt.display = 0;
           p.model = minFunc(funObjL2,w_init,opt);
       
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
         % return the parameters 
         
         m = {obj.params.model(2:end)'}; % ignore bias term
         desc = {'unknown'};
                    
       end
       
    end
end 
