classdef lrgauss_adf < classifier
%Posterior for a Gaussian prior and a logistic likelihood
%
%   Options:
%   'prior' : the prior Gaussian distribution in the format
%              prior = struct('mean',mean_vector,'cov',cov_matix)
%           : if none uses as prior the standard Gaussian    
%
%   SEE ALSO:
%   lrgauss_lap
%
%   Copyright (c) 2008, Adriana Birlutiu
%
%   $Log: lrgauss_adf.m,v $
%

    properties

      data;
      targets;
      
      model; % the weight vector
      prior;
        
    end

    methods
       function obj = lrgauss_adf(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
            
           if iscell(data), error('GP does not take multiple datasets as input'); end
            
           if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
                      
           %if obj.nclasses ~= 2, error('lrgauss_adf only accepts binary class problems'); end
                   
           targets = design(:,1);
           targets(design == 1) = 1;
           targets(design == 2) = -1;

           if isempty(obj.prior)
               prior = priorstandard(size(data,2),1);
           else
               prior = obj.prior;
           end
           
           % training mode
           prior = logisticgauss_adf(prior, data, targets);
           obj.model = [ prior.mean' 0; -prior.mean' 0 ];
           
           obj.data = data;
           obj.targets = targets;
       end
       
        function post = test(obj,data)
           
           if iscell(data)
               post = cell(1,length(data));
               for j=1:length(data)
                   post{j} = slr_classify([data{j} ones(size(data{j},1),1)], obj.model{j});
               end
           else       
               post = slr_classify([data ones(size(data,1),1)], obj.model);
           end
       end          
      
    end
end 
