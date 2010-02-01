classdef static_classifier < classifier
%STATIC_CLASSIFIER class to represent classifiers that use a Bayesian
%network structure as their underlying basis.
%
%   The variable of interest is always assumed to be the first variable in
%   each time slice.
%
%   NOTE:
%       this is an abstract class which requires another classifier
%       that specifies the structure of the associated graphical model
%
%   SEE ALSO:
%       
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: static_classifier.m,v $
%

    properties
        ie;
        bn;         % save for info
        engine;     % this engine is used if specified
        numvar;     % number of variables
        nclasses;   % number of classes
    end

    methods
      
       function obj = static_classifier(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = estimate(obj,data,design)
            
            p.nclasses = design.nunique; 

            data = data.X;
            design = design.X;
            
            if isempty(obj.numvar), obj.numvar = size(data,2)+1; end
            
            % construct factors (custom code)
            factors = obj.construct_factors();

            % construct BN
            bn = bayesnet(factors);
                        
            % learn parameters using the standard learner
            bn = bn.learn_parameters([design(:,1) data]);
                           
            % create inference engine
            if isempty(obj.engine)
                p.ie = hugin_ie(bn);
            else
                p.ie = obj.engine(bn);
            end
                                      
            p.bn = bn;
            
       end
       
       function post = map(obj,data)       
           
           data = data.X;
           
           post = zeros([size(data,1) obj.params.nclasses]);
           
           for j=1:size(post,1)

               % add evidence to the inference engine
               obj.params.ie.enter_evidence([nan data(j,:)]);

               % compute marginal for first variable
               m = normalize(obj.params.ie.marginalize(1));

               post(j,:) = m.p';
           end

           post = dataset(post);
           
       end

        function factors = construct_factors(obj)
        % This function should be overloaded. It is currently used 
        % to return the factors of a prespecified DBN.
        
        if isempty(obj.params.bn)
          error('unspecified BN');
        end
        
        factors = obj.params.bn.factors;
        
       end
       
    end
end 
