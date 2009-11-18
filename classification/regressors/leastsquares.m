classdef leastsquares < regressor
%LEASTSQUARES regressor
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: leastsquares.m,v $
%

    properties
        
      model     
      
    end

    methods
       function obj = leastsquares(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function obj = train(obj,data,design)
        
         % add bias term
         data = [data ones(size(data,1),1)];         
         obj.model = data\design;
         
       end
       
       function post = test(obj,data)       
       
         % add bias term
         data = [data ones(size(data,1),1)];
       
         post = data * obj.model;

       end
       
       function m = getmodel(obj,label,dims)
         % return the parameters wrt a class label in some shape 
                    
         m = full(obj.model(1:(end-1))); % ignore bias term

         if nargin == 3 && numel(m) == prod(dims)
           m = reshape(m,dims);
         end
         
       end

    end
end 
