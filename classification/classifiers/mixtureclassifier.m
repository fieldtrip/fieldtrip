classdef mixtureclassifier < static_classifier
%MIXTURECLASSIFIER class 
%
%   assumes as structure f <- m -> c where the f's are feature; m is a
%   mixture and c is the class variable. The number of mixture components
%   is specified by mixture. If unspecified the number of mixture
%   components equals the number of classes
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: mixtureclassifier.m,v $
%

    properties
        mixture; % number of mixture components
    end

    methods
       function obj = mixtureclassifier(varargin)
                  
           obj = obj@classifier(varargin{:});
           
       end
       function obj = train(obj,data,design)

         [data,design] = obj.check_input(data,design);
         
           if isempty(obj.mixture), obj.mixture = max(design(:,1)); end

           % data must accommodate hidden variable         
           obj = obj.train@static_classifier([data nan(size(data,1),1)],design);
            
       end
       
       function post = test(obj,data)

         data = obj.check_input(data);
         
         post = obj.test@static_classifier([data nan(size(data,1),1)]);
       
       end

       function factors = construct_factors(obj)
                 
           % definition of factors

           factors = cell(1,obj.numvar);

           % class variable
           factors{1} = multinomial_cpd(1,obj.numvar,ones([obj.nclasses obj.mixture]));
           
           % feature variables
           for f=2:(obj.numvar-1)
               factors{f} = gaussian_cpd(f,[],obj.numvar,randn([obj.mixture 1]),cell([obj.mixture 1]),rand([obj.mixture 1]));
           end

           % hidden variable
           factors{obj.numvar} = multinomial_cpd(obj.numvar,[], 1/obj.mixture - 1e-2 + 2e-2*rand([obj.mixture 1]));               
       end            
       
    end

end
