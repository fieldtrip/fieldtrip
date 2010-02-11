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
       
       function p = estimate(obj,X,Y)
                 
         nclasses = obj.nunique(Y);
         
         w = smlr_learning(Y, [X ones(size(X,1),1)], size(X,2)+1, ...
           'wdisplay', 'off', 'wmaxiter', 50', 'nlearn', 300, 'nstep', 100,...
           'amax', 1e8, 'isplot', 0, 'gamma0', 0);

         p.model = reshape(w, [size(X,2)+1, nclasses]);

       end

       function Y = map(obj,X)

         [tmp, Y] = calc_label([X ones(size(X,1),1)], obj.params.model);
         
       end              
       
       function [m,desc] = getmodel(obj)
         % return the parameters wrt a class label in some shape 
         
         m = mat2cell(obj.params.model(:,1:(end-1))',ones(1,size(obj.params.model(:,1:(end-1))',1)),size(obj.params.model(:,1:(end-1))',2));
         desc = {'unknown'};
           
       end
       
    end
end 
