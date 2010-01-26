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

        % model
        data;
        design;
        
        % options
        disfun = 'diagLinear';
        
        % diagnostics
        coeff;
        
    end

    methods
       function obj = da(varargin)
       
           % check availability
           if ~license('test','statistics_toolbox')
               error('requires Matlab statistics toolbox');
           end
           
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
            % simply stores input data and design
            
            obj.data = data;
            obj.design = design;
                       
       end
       function post = test(obj,data)       

           [class,err,post,logp,obj.coeff] = classify(data.X,obj.data.X,obj.design.X,obj.disfun);                    
       
           post = dataset(post);
       end

    end
end 
