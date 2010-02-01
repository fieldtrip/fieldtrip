classdef nearestneighbour < classifier
%nearest neighbour classifier
%
% SEE ALSO:
% knnclassify
%
% REQUIRES:
% bioinformatics toolbox
%
%   Copyright (c) 2008, Marcel van Gerven
%

    properties
  
      k=1; % number of neighbours
      distance = 'euclidean'; % distance function
      rule = 'nearest'; % classification rule
        
    end

    methods
      
       function obj = nearestneighbour(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function p = map(obj,data,design)
         % simply store training data
         
         p.data = data;
         p.design = design;
         
       end
       
       function post = estimate(obj,data)
         
         l = knnclassify(data.X,obj.params.data.X,obj.params.design.X,obj.k,obj.distance,obj.rule);
         
         post = zeros(size(data,1),obj.params.design.nunique);
         for j=1:data.nsamples
           post(j,l(j)) = 1;
         end
         
         post = dataset(post);
         
       end
       
    end
end
