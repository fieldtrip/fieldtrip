classdef ridge < regressor
%RIDGE regressor
%
% NOTE: 
% using elasticnet with lambda = 0 and nu ~= 0 can be faster
%
%   Copyright (c) 2009, Marcel van Gerven


    properties
        
      lambda = 1;
      
    end

    methods
       function obj = ridge(varargin)
          
         obj = obj@regressor(varargin{:});
         
       end
       function p = estimate(obj,X,Y)
        
         % add bias term
         X = [X ones(size(X,1),1)];
          
         lambdas = obj.lambda*ones(size(X,2),1); % Penalize the absolute value of each element by the same amount

         R = chol(X'*X + diag(lambdas));
         
         p.model = R\(R'\(X'*Y));
         
       end
       
       function Y = map(obj,X)       
       
         Y = [X ones(size(X,1),1)] * obj.params.model;

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters
                    
         m = {full(obj.params.model(1:(end-1)))}; % ignore bias term
         desc = {'unknown'};
         
       end

    end
end 
