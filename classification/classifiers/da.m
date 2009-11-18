classdef da < classifier
%DA linear discriminant analysis
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
            
            [data,design] = obj.check_input(data,design);
            
            if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
            
            obj.data = data;
            obj.design = design(:,1);
                       
       end
       function post = test(obj,data)       

         data = obj.check_input(data);
            
           % deal with empty data
           if size(data,2) == 0
               
               % random assignment
               post = rand([size(data,1) obj.nclasses]);               
               post = double(post == repmat(max(post,[],2),[1 obj.nclasses]));
               return
           end

           [class,err,post,logp,obj.coeff] = classify(data,obj.data,obj.design(:,1),obj.disfun);                    
       end

    end
end 
