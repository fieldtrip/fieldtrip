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
       
       function p = estimate(obj,X,Y)
                 
         if obj.nunique(Y) ~= 2, error('SLRLAP expects binary class labels'); end

         p.model = slr_learning_var2(Y-1, [X ones(size(X,1),1)],...
           'nlearn', 300, 'nstep', 100, 'amax', 1e8, 'invhessian', obj.invhessian);


       end

       function Y = map(obj,X)

         [tmp, Y] = calc_label([X ones(size(X,1),1)], [zeros(size(X,2)+1,1) obj.params.model]);
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.params.model(1:(end-1))'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
