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
            
           if iscell(data), error('GP does not take multiple datasets as input'); end
            
           
           if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end

           
         %  if obj.nclasses ~= 2, error('lrgauss_lap only accepts binary class problems'); end
           
           targets = design(:,1);
           targets(design == 1) = -1;
           targets(design == 2) = 1;

           % training mode
           x = logisticgauss_lap(obj.prior, data, targets);
           obj.model = [ [x' 0]; [-x' 0] ];
             
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
