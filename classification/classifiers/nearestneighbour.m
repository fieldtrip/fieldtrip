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
       
       function p = estimate(obj,X,Y)
         % simply store training data
         
         p.data = X;
         p.design = Y;
         
       end
       
       function Y = map(obj,X)
         
         l = knnclassify(X,obj.params.data,obj.params.design,obj.k,obj.distance,obj.rule);
         
         Y = zeros(size(X,1),obj.nunique(obj.params.design));
         for j=1:size(X,1)
           Y(j,l(j)) = 1;
         end
         
       end
       
    end
end
