classdef da < classifier
% DA linear discriminant analysis
%
%   Options:
%   'disfun' : type of discriminant function used
%
%   SEE ALSO:
%   classify.m
%
%   REQUIRES:
%   Matlab(R) Statistics Toolbox
%
%   Copyright (c) 2008, Marcel van Gerven


    properties

        % options
        disfun = 'diagLinear';        
        
    end

    methods
      
       function obj = da(varargin)
       
           % check availability
           if ~license('test','statistics_toolbox')
               error('requires Matlab statistics toolbox');
           end
           
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = estimate(obj,X,Y)
            % simply stores input data and design
            
            p.data = X;
            p.design = Y;
                       
       end
       
       function post = map(obj,X)       

           [class,err,post] = classify(X,obj.params.data,obj.params.design,obj.disfun);                    
       
       end

    end
end 
