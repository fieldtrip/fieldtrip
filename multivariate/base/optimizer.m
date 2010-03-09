classdef optimizer < predictor
%OPTIMIZER can be used to optimize certain parameters of a method
%
%   Usage:
%
%    optimizer('mvmethod',l1lr,'validator',crossvalidator(),'vars','lambda','vals',[1 10 100],'metric','accuracy')
%
%     will act as if it is l1lr method but will optimize the variable
%     lambda in the range [1 10 100] using the specified validator and metric.
%
%   parameters:
%    configs               % variables configurations
%    results               % configuration results
%    optimum               % optimal configuration
%
%   Remarks:
%     this replaces the obsolete optimize.m function
%     mvmethods are supposed to handle optimization efficiently by starting
%     at previous solutions specified in the .params field; l1lr provides
%     an example
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

      mvmethod              % the method for which to optimize
      
      validator             % the used validator
      
      vars                  % the variables to optimize
      vals                  % the values used
      metric = 'accuracy';  % the evaluation metric  
                  
    end

    methods
      
       function obj = optimizer(varargin)
         
         % parse options
         for i=1:2:length(varargin)
           obj.(varargin{i}) = varargin{i+1};
         end
         
         if isempty(obj.validator)
           obj.validator = crossvalidator('verbose',true);
         end
         assert(isa(obj.validator,'validator'));
         
         assert(~isempty(obj.mvmethod));
         assert(~isempty(obj.vars));
         assert(~isempty(obj.vals));
         
         if ~iscell(obj.vars)
           obj.vars = { obj.vars };
           obj.vals = { obj.vals };
         end
         
         obj.validator.procedure = mva({obj.mvmethod});
         
       end
       
       function p = estimate(obj,data,design)
         
         if obj.verbose
           fprintf('optimizing variables');
           for c=1:length(obj.vars)
             fprintf(' %s',obj.vars{c});
           end
           fprintf('\n');
         end
         
         % the optimizer iterates over all values for the specified
         % variables and returns the best result

         nv = length(obj.vars);

         maxresult = -Inf;

         vld = obj.validator; % the procedure used by the validator
      
         % iterate over configurations
         p.configs = cartprod(obj.vals{:});
         for i=1:size(p.configs,1)

           % set parameters
           if iscell(vld.procedure)
           % happens when the vld has already been validated
             for c=1:numel(vld.procedure)
               for j=1:nv
                vld.procedure{c}.mvmethods{1}.(obj.vars{j}) = p.configs(i,j);
               end
             end
           else
             for j=1:nv
               vld.procedure.mvmethods{1}.(obj.vars{j}) = p.configs(i,j);
             end
           end
           
           if obj.verbose
             fprintf('evaluating configuration ');
             for j=1:size(p.configs,2)
               fprintf(' %d',p.configs(i,j));
             end
             fprintf('\n');
           end

           vld = vld.validate(data,design);
           p.results(i) = vld.evaluate('metric',obj.metric);

           if obj.verbose
             fprintf('%s : %f\n',obj.metric,p.results(i));
           end
           
           if p.results(i) > maxresult

             maxresult = p.results(i);
             p.optimum = i; % optimal configuration
           end

         end
         
         if obj.verbose
           fprintf('optimum configuration ');
           for j=1:size(p.configs,2)
             fprintf(' %d',p.configs(p.optimum,j));
           end
           fprintf(' : %f\n',p.results(p.optimum));
         end
         
         % we now know the optimal configuration. Now we retrieve the
         % method of interest with this configuration and retrain
      
         p.mvmethod = obj.mvmethod;
         for j=1:nv
           p.mvmethod.(obj.vars{j}) = p.configs(p.optimum,j);
         end
         p.mvmethod = p.mvmethod.train(data,design);
         
       end
       
       function post = map(obj,data)       
           
         post = obj.params.mvmethod.map(data);
         
       end
       
       function p = predict(obj,data)
        % relevant in case the optimizer acts as a predictor
        
        if isa(obj.method,'predictor')
          p = obj.mvmethod.predict(data);
        else
          if obj.verbose, fprintf('%s cannot predict\n',class(obj.mvmethod)); end
          p = nan(size(data,1),1);
        end
         
       end
       
       function metric = evaluate(obj,X,Y,varargin)
         
         metric = obj.mvmethod.evaluate(X,Y,varargin{:});
         
       end
       
       function  p = significance(obj,X,Y,varargin)
         
         p = obj.mvmethod.significance(X,Y,varargin{:});
         
       end
       
       function [m,desc] = getmodel(obj)
         % call the enclosing method
         
         [m,desc] = obj.params.mvmethod.getmodel();
           
       end

    end
end 
