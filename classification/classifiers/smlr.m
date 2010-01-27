classdef smlr < classifier
%SMLR wrapper to sparse multinomial logistic regression
%
% http://www.cns.atr.jp/~oyamashi/SLR_WEB.html
%
% Copyright (c) 2010, Marcel van Gerven


    properties

        model; 
        
    end

    methods
       
      function obj = smlr(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function obj = train(obj,data,design)
                 
         nclasses = design.nunique;
         
         w = smlr_learning(design.X, [data.X ones(data.nsamples,1)], data.nfeatures+1, ...
           'wdisplay', 'off', 'wmaxiter', 50', 'nlearn', 300, 'nstep', 100,...
           'amax', 1e8, 'isplot', 0, 'gamma0', 0);

         obj.model = reshape(w, [data.nfeatures+1, nclasses]);

       end

       function post = test(obj,data)

         [tmp, post] = calc_label([data.X ones(data.nsamples,1)], obj.model);
         
         post = dataset(post);
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
           m = {obj.model(2:end)'}; % ignore bias term
           desc = {'unknown'};
           
       end
       
    end
end 
