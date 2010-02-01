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
%
%   $Log: da.m,v $
%

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
       
       function p = estimate(obj,data,design)
            % simply stores input data and design
            
            p.data = data;
            p.design = design;
                       
       end
       
       function post = map(obj,data)       

           [class,err,post] = classify(data.X,obj.params.data.X,obj.params.design.X,obj.disfun);                    
       
           post = dataset(post);
       end

    end
end 
