classdef lrgauss_lap < classifier
%Posterior for a Gaussian prior and a logistic likelihood
%
%   Options:
%   'prior' : the prior Gaussian distribution in the format
%              prior = struct('mean',mean_vector,'cov',cov_matix)
%           : if none maximizes only the likelihood    
%
%   SEE ALSO:
%   lrgauss_adf
%
%   Copyright (c) 2008, Adriana Birlutiu
%
%   $Log: lrgauss_lap.m,v $
%

    properties

      data;
      model; % the weight vector
      prior;
        
    end

    methods
       function obj = lrgauss_lap(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
          
         if design.nunique ~= 2, error('LRGAUSS_LAP only accepts binary class problems'); end
                  
         targets = design.X(:,1);
         targets(targets == 1) = -1;
         targets(targets == 2) = 1;
         
         % training mode
         x = logisticgauss_lap(obj.prior, data.X, targets);
         obj.model = [ [x' 0]; [-x' 0] ];
         
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
