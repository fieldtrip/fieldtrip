classdef gpregressor < regressor
%GP gaussian process regressor
%
%   Options:
%   'optimize' : if true tries to optimize hyperparameters
%
%   SEE ALSO:
%   gpml-matlab
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

        data;
        targets;
        
        optimize = true;
        loghyper; % log likelihood of the hyperparameters
        covfunc = {'covSum', {'covSEard','covNoise'}}; % covariance function
        offset; % offset from zero for the targets

    end

    methods
      
       function obj = gpregressor(varargin)
                  
           obj = obj@regressor(varargin{:});
           
       end
       
       function p = estimate(obj,data,design)
            
           if ~exist('gpml-matlab','dir')
               error('this code requires an external toolbox: http://www.gaussianprocess.org/gpml/code/matlab/doc/');
           end
           
           targets = design.X(:,1);
           
           % center targets
           p.offset = mean(targets);
           targets = targets - p.offset;
           
           p.data = data.X;
           p.targets = targets;

          if obj.optimize % optimize hyperparameters
             p.loghyper = minimize([zeros(1,data.nfeatures) 0 log(sqrt(0.1))]', 'gpr', -100, obj.covfunc, data.X, targets);
           end
                       
       end
       
       function post = map(obj,data)       
           % returns mean and variance
       
           [avg,variance] = gpr(obj.params.loghyper, obj.covfunc, obj.params.data, obj.params.targets, data.X);
           avg = avg + obj.params.offset;  % add back offset to get true prediction

           post = dataset([avg variance]);
           
       end

    end
end 
