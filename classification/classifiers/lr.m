classdef lr < classifier
% logistic regression using Mark Schmidt's L1General package
%
%   Copyright (c) 2010, Marcel van Gerven


    properties

        model; 
        
    end

    methods
       function obj = lr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
                 
         if design.nunique > 2, error('LR expects binary class labels'); end
         
         if obj.verbose
           fprintf('\nComputing Logistic Regression Coefficients...\n');
         end

         design = design.collapse();
         
         % foolproof
         if all(design(:,1) == 1)
           obj.model = -inf;
         elseif all(design(:,1) == 2)
           obj.model = inf;
         else
           
           X = [ones(data.nsamples,1) data.collapse()];
           
           design = 2*design - 3; % convert to -1 / + 1
           
           funObj = @(w)LogisticLoss(w,X,design);
           w_init = zeros(data.nfeatures+1,1);
           
           mfOptions.Method = 'lbfgs';
           mfOptions.display = 0;
           obj.model = minFunc(funObj,w_init,mfOptions);
         
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
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
         m = {obj.model(2:end)'}; % ignore bias term
         desc = {'unknown'};
           
       end
       
    end
end 
