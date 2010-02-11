classdef linreg < regressor
%LINREG linear regression method class
%
%   Copyright (c) 2008, Marcel van Gerven

    
    methods
        
        function obj = linreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function p = estimate(obj,X,Y)
            
            p.model = regress(Y,[X ones(size(X,1),1)]);            
        end
        
        function Y = map(obj,X)           
            Y = [X ones(size(X,1),1)] * obj.params.model;            
        end
 
    end
end 
