classdef leastsquares < regressor
%LEASTSQUARES regressor
%
%   Copyright (c) 2009, Marcel van Gerven


    methods
       
      function obj = leastsquares(varargin)
          
         obj = obj@regressor(varargin{:});
         
      end
      
      function p = estimate(obj,X,Y)
         
         p.model = [X ones(size(X,1),1)]\Y;
         
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
