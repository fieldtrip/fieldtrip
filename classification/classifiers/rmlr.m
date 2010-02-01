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
       
       function p = estimate(obj,data,design)
                 
         nclasses = design.nunique;
         
         w = rmlr_learning(design.X, [data.X ones(data.nsamples,1)], data.nfeatures+1, ...
           'wdisplay', 'off', 'wmaxiter', 50', 'nlearn', 300, 'nstep', 100,...
           'amax', 1e8, 'gamma0', 0);
     
         p.model = reshape(w, [data.nfeatures+1, nclasses]);

       end

       function post = map(obj,data)

         [tmp, post] = calc_label([data.X ones(data.nsamples,1)], obj.params.model);
         
         post = dataset(post);
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.model(2:end)'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
