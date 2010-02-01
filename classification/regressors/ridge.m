classdef ridge < regressor
%RIDGE regressor
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: ridge.m,v $
%

    properties
        
      lambda = 1;
      
    end

    methods
       function obj = ridge(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function p = estimate(obj,data,design)
        
         % add bias term
         data = [data.X ones(data.nsamples,1)];
         design = design.X;
         
         lambdas = obj.lambda*ones(size(data,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(data'*data + diag(lambdas));
         
         p.model = R\(R'\(data'*design));
         
       end
       
       function post = map(obj,data)       
       
         post = dataset([data.X ones(data.nsamples,1)] * obj.params.model);

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters
                    
         m = {full(obj.params.model(1:(end-1)))}; % ignore bias term
         desc = {'unknown'};
         
       end

    end
end 
