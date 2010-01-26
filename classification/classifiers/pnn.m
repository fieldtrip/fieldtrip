classdef pnn < classifier
%PNN probabilistic neural network classifier
%
%   Options
%   'spread' : newpnn parameter that will be optimized if it is a vector
%
%   SEE ALSO:
%   newpnn.m
%
%   NOTE:
%   The FieldTrip replacement of dist.m does not work in this context!
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: pnn.m,v $
%

    properties
        net; % the neural network
        spread = 0.1; % modifies nn behaviour; if a vector, will be optimized
    end

    methods
      function obj = pnn(varargin)
        
        obj = obj@classifier(varargin{:});
        
        % check availability
        if ~license('test','neural_network_toolbox')
          error('requires Matlab neural network toolbox');
        end
        
      end
      
      function obj = train(obj,data,design)
        
        obj.net = newpnn(data.X',design.X(:,1)',obj.spread);
        
      end
      
      function post = test(obj,data)
        
        post = dataset(sim(obj.net,data.X')');
      end
      
    end
end
