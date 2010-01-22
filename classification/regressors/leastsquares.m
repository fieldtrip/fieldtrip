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
         
         obj.model = [data.collapse() ones(data.nsamples,1)]\design.collapse();
         
       end
       
       function post = test(obj,data)       
           
         post = dataset([data.collapse() ones(data.nsamples,1)]) * obj.model;

       end
       
       function m = getmodel(obj)
         % return the parameters 
                    
         m = full(obj.model(1:(end-1))); % ignore bias term
   
       end

    end
end 
