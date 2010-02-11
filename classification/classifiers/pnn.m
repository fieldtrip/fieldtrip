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
      
      function p = estimate(obj,X,Y)
        
        p.net = newpnn(X',Y(:,1)',obj.spread);
        
      end
      
      function Y = map(obj,X)
        
        Y = sim(obj.params.net,X')';
      end
      
    end
end
