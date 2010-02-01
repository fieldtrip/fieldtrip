classdef linreg < regressor
%LINREG linear regression method class
%
%   Copyright (c) 2008, Marcel van Gerven

    
    methods
        
        function obj = linreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function p = estimate(obj,data,design)
            
            p.model = regress(design.X,[data.X ones(data.nsamples,1)]);            
        end
        
        function res = map(obj,data)           
            res = dataset([data.X ones(data.nsamples,1)] * obj.params.model);            
        end
 
    end
end 
