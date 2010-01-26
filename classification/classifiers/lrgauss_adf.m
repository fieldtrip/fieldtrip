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
        
         if design.nunique ~= 2, error('LRGAUSS_ADF only accepts binary class problems'); end

         targets = design.X(:,1);
         targets(targets == 1) = -1;
         targets(targets == 2) = 1;
         
         if isempty(obj.prior)
           prior = priorstandard(data.nfeatures,1);
         else
           prior = obj.prior;
         end
         
         % training mode
         prior = logisticgauss_adf(prior, data.X, targets);
         obj.model = [  -prior.mean' 0; prior.mean' 0];
         
         obj.data = data.X;
         obj.targets = targets;
       end
       
       function post = test(obj,data)
         
         if iscell(data)
           post = cell(1,length(data));
           for j=1:length(data)
             post{j} = slr_classify([data{j}.X ones(data{j}.nsamples,1)], obj.model{j});
           end
         else
           post = slr_classify([data.X ones(data.nsamples,1)], obj.model);
         end
         
         post = dataset(post);
         
       end
       
    end
end
