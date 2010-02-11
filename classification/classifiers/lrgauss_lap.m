classdef lrgauss_lap < classifier
%Posterior for a Gaussian prior and a logistic likelihood
%
%   Options:
%   'prior' : the prior Gaussian distribution in the format
%              prior = struct('mean',mean_vector,'cov',cov_matix)
%
%   SEE ALSO:
%   lrgauss_adf
%
%   Copyright (c) 2010, Adriana Birlutiu


    properties

      prior;
        
    end

    methods
       function obj = lrgauss_lap(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function p = estimate(obj,X,Y)
          
         if obj.nunique(Y) ~= 2, error('LRGAUSS_LAP only accepts binary class problems'); end
                  
         targets = Y(:,1);
         targets(targets == 1) = -1;
         targets(targets == 2) = 1;
         
         if isempty(obj.prior)
           pr = priorstandard(size(X,2)+1,1);
         else
           pr = obj.prior;
         end
         
         % training mode
         x = logisticgauss_lap(pr, [X ones(size(X,1),1)], targets);
         p.model = [ x'; -x' ];
         
       end
       
       function Y = map(obj,X)
        
         Y = slr_classify([X ones(size(X,1),1)], obj.params.model);
         
       end
       
    end
end
