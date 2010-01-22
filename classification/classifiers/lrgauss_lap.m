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
         
         data = data.collapse();
         design = design.collapse();
         
         targets = design(:,1);
         targets(design == 1) = -1;
         targets(design == 2) = 1;
         
         % training mode
         x = logisticgauss_lap(obj.prior, data, targets);
         obj.model = [ [x' 0]; [-x' 0] ];
         
       end
       
       function post = test(obj,data)
         
         data = data.collapse();
         
         if iscell(data)
           post = cell(1,length(data));
           for j=1:length(data)
             post{j} = slr_classify([data{j} ones(size(data{j},1),1)], obj.model{j});
           end
         else
           post = slr_classify([data ones(size(data,1),1)], obj.model);
         end
         
         post = dataset(post);
         
       end
       
    end
end
