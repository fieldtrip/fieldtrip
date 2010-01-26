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
         
         targets = zeros(data.nsamples,obj.nclasses);
         for j=1:data.nsamples
           targets(j,design.X(j,1)) = 1;
         end
         
         obj.net = knn(data.nfeatures, obj.nclasses, obj.k, data.X, targets);
         
       end
       
       function post = test(obj,data)
         
        [y, l] = knnfwd(obj.net, data.X);
         
         post = zeros(size(data,1),obj.nclasses);
         for j=1:data.nsamples
           post(j,l(j)) = 1;
         end
         
         post = dataset(post);
         
       end
       
    end
end
