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
%
%   $Log: gpregressor.m,v $
%

    properties

        data;
        targets;
        
        optimize = true;
        loghyper; % log likelihood of the hyperparameters
        covfunc; % covariance function
        offset; % offset from zero for the targets

    end

    methods
       function obj = gpregressor(varargin)
                  
           obj = obj@regressor(varargin{:});
                      
       end
       function obj = train(obj,data,design)
            
           if iscell(data), error('GPREGRESSOR does not take multiple datasets as input'); end
            
           if ~exist('gpml-matlab','dir')
               error('this code requires an external toolbox: http://www.gaussianprocess.org/gpml/code/matlab/doc/');
           end
           
           targets = design(:,1);
           
           % center targets
           obj.offset = mean(targets);
           targets = targets - obj.offset;
           
           obj.data = data;
           obj.targets = targets;

           % specify standard covariance function
           obj.covfunc = {'covSum', {'covSEard','covNoise'}};

           if obj.optimize % optimize hyperparameters
             obj.loghyper = minimize([zeros(1,size(data,2)) 0 log(sqrt(0.1))]', 'gpr', -100, obj.covfunc, data, targets);
           end
                       
       end
       function post = test(obj,data)       
           % returns mean and variance
           
           if iscell(data), error('GPREGRESSOR does not take multiple datasets as input'); end

           [avg,variance] = gpr(obj.loghyper, obj.covfunc, obj.data, obj.targets, data);
           avg = avg + obj.offset;  % add back offset to get true prediction

           post = [avg variance];
           
       end

    end
end 
