classdef rmlr < classifier
%RMLR wrapper to regularized multinomial logistic regression
%
% http://www.cns.atr.jp/~oyamashi/SLR_WEB.html
%
% Copyright (c) 2010, Marcel van Gerven


    properties

        model; 
        
    end

    methods
       
      function obj = rmlr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = estimate(obj,X,Y)
                 
         nclasses = obj.nunique(Y);
         
         w = rmlr_learning(Y, [X ones(size(X,1),1)], size(X,2)+1, ...
           'wdisplay', 'off', 'wmaxiter', 50', 'nlearn', 300, 'nstep', 100,...
           'amax', 1e8, 'gamma0', 0);
     
         p.model = reshape(w, [size(X,2)+1, nclasses]);

       end

       function post = map(obj,X)

         [tmp, post] = calc_label([X ones(size(X,1),1)], obj.params.model);
          
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.model(2:end)'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
