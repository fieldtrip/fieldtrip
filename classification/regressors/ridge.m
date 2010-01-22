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
         data = [data.collapse() ones(data.nsamples,1)];
         design = design.collapse();
         
         lambdas = obj.lambda*ones(size(data,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(data'*data + diag(lambdas));
         
         obj.model = R\(R'\(data'*design));
         
       end
       
       function post = test(obj,data)       
       
         post = dataset([data.collapse() ones(data.nsamples,1)] * obj.model);

       end
       
       function m = getmodel(obj)
         % return the parameters
                    
         m = full(obj.model(1:(end-1))); % ignore bias term
         
       end

    end
end 
