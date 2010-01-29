classdef slrvar < classifier
%SLRVAR wrapper to sparse logistic regression with a variational approximation
%
% http://www.cns.atr.jp/~oyamashi/SLR_WEB.html
%
% Copyright (c) 2010, Marcel van Gerven


    properties

        model; 
        
        invhessian = 0; % 1 is faster but too expensive for many features
        
    end

    methods
       
      function obj = slrvar(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function obj = train(obj,data,design)
                 
         if design.nunique ~= 2, error('SLRLAP expects binary class labels'); end

         design = design.X-1; % zero based
         
         obj.model = slr_learning_var2(design, [data.X ones(data.nsamples,1)],...
           'nlearn', 300, 'nstep', 100, 'amax', 1e8, 'invhessian', obj.invhessian);


       end

       function post = test(obj,data)

         [tmp, post] = calc_label([data.X ones(data.nsamples,1)], [zeros(data.nfeatures+1,1) obj.model]);
         
         post = dataset(post);
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.model(2:end)'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
