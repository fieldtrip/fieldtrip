classdef ridge < regressor
%RIDGE regressor
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: ridge.m,v $
%

    properties
        
      model     
      lambda = 1;
      
    end

    methods
       function obj = ridge(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function obj = train(obj,data,design)
        
         % add bias term
         data = [data ones(size(data,1),1)];
         
         lambdas = obj.lambda*ones(size(data,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(data'*data + diag(lambdas));
         
         obj.model = R\(R'\(data'*design));
         
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
