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

      data
      design
      
      k=1; % number of neighbours
      distance = 'euclidean'; % distance function
      rule = 'nearest'; % classification rule
        
    end

    methods
      
       function obj = nearestneighbour(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       
       function obj = train(obj,data,design)
         % simply store training data
         
         obj.data = data;
         obj.design = design;
         
       end
       
       function post = test(obj,data)
         
         l = knnclassify(data.X,obj.data.X,obj.design.X,obj.k,obj.distance,obj.rule);
         
         post = zeros(size(data,1),obj.design.nunique);
         for j=1:data.nsamples
           post(j,l(j)) = 1;
         end
         
         post = dataset(post);
         
       end
       
    end
end
