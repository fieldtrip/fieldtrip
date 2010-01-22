classdef nearestneighbour < classifier
%NEAREST NEIGHBOUR classifier
%
%  Requires Netlab
%
%   SEE ALSO:
%   knn
%   knnfwd
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: nearestneighbour.m,v $
%

    properties

        k=1; % number of neighbours
        net; % the knn object        
        
        nclasses;
        
    end

    methods
       function obj = nearestneighbour(varargin)
                  
           obj = obj@classifier(varargin{:});
                      
       end
       function obj = train(obj,data,design)
         % create knn object
         
         obj.nclasses = design.nunique;
         
         data = data.collapse();
         design = design.collapse();
         
         targets = zeros(size(data,1),obj.nclasses);
         for j=1:size(data,1)
           targets(j,design(j,1)) = 1;
         end
         
         obj.net = knn(size(data,2), obj.nclasses, obj.k, data, targets);
         
       end
       function post = test(obj,data)
         
         data = data.collapse();
         
         [y, l] = knnfwd(obj.net, data);
         
         post = zeros(size(data,1),obj.nclasses);
         for j=1:size(data,1)
           post(j,l(j)) = 1;
         end
         
         post = dataset(post);
         
       end
       
    end
end
