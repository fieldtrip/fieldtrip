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
         
         obj.model = [data.X ones(data.nsamples,1)]\design.X;
         
       end
       
       function post = test(obj,data)       
           
         post = dataset([data.X ones(data.nsamples,1)]) * obj.model;

       end
       
       function [m,desc] = getmodel(obj)
         % return the parameters 
                    
         m = {full(obj.model(1:(end-1)))}; % ignore bias term
         desc = {'unknown'};
         
       end

    end
end 
