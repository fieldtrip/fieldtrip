classdef linreg < regressor
%LINREG linear regression method class
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: linreg.m,v $
%
    properties        
        
        model % regression coefficients        
    end
    
    methods
        
        function obj = linreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function obj = train(obj,data,design)
            
            obj.model = regress(design.collapse(),[data.collapse() ones(data.nsamples,1)]);            
        end
        
        function res = test(obj,data)           
            res = [data.collapse() ones(data.nsamples,1)] * obj.model;            
        end
 
    end
end 
