classdef mixtureclassifier < static_classifier
%MIXTURECLASSIFIER class 
%
%   assumes as structure f <- m -> c where the f's are feature; m is a
%   mixture and c is the class variable. The number of mixture components
%   is specified by mixture. If unspecified the number of mixture
%   components equals the number of classes
%
%   Copyright (c) 2008, Marcel van Gerven


    properties
        mixture; % number of mixture components
    end

    methods
       function obj = mixtureclassifier(varargin)
                  
           obj = obj@static_classifier(varargin{:});
           
       end
       
       function p = estimate(obj,X,Y)
  
         if isempty(obj.mixture), obj.mixture = obj.nunique(Y); end
         
         % data must accommodate hidden variable
         p = obj.estimate@static_classifier([X nan(size(X,1),1)],Y);
         
       end
       
       function Y = map(obj,X)

         Y = obj.map@static_classifier([X nan(size(X,1),1)]);
       
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
