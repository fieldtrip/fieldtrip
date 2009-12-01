classdef optimizer < predictor
%OPTIMIZER can be used to optimize certain parameters of a method
%
%   Usage:
%   optimizer('validator',crossvalidator('procedure',
%     {kernelmethod()}),'variables','C','values',1:10,'metric','accuracy')
%
%     will act as if it is a kernelmethod but will optimize the variable C
%     in the range 1:10 using the specified validator and metric.
%
%   Remarks:
%     this replaces the obsolete optimize.m function
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: optimizer.m,v $
%

    properties

      validator             % the used validator
      variables             % the variables to optimize
      values                % the values used
      metric = 'accuracy';  % the evaluation metric  
      methodidx = 1;        % index of the method in the validator procedure that is to be optimized (default = 1)
      
      configs               % variables configurations
      results               % configuration results
      optimum               % optimal configuration
      
      method                % the method of interest
      
    end

    methods
       function obj = optimizer(varargin)
       
          % parse options
          for i=1:2:length(varargin)
            obj.(varargin{i}) = varargin{i+1};
          end

         assert(isa(obj.validator,'validator'));

         assert(~isempty(obj.variables));
         assert(~isempty(obj.values));

         if ~iscell(obj.variables)
           obj.variables = { obj.variables };
           obj.values = { obj.values };
         end
                               
       end
       function obj = train(obj,data,design)
         
         if obj.verbose
           fprintf('optimizing\n');
         end
         
         % the optimizer iterates over all values for the specified
         % variables and returns the best result

         nv = length(obj.variables);

         maxresult = -Inf;

         % iterate over configurations
         obj.configs = cartprod(obj.values{:});
         for i=1:size(obj.configs,1)

           vld = obj.validator; % the procedure used by the validator
           
           % set parameters
           for j=1:nv
             vld.procedure.clfmethods{obj.methodidx}.(obj.variables{j}) = obj.configs(i,j);
           end

           if obj.verbose
             fprintf('evaluating configuration ');
             fprintf('%d ',obj.configs(i,:));
           end

           v = vld.validate(data,design);
           obj.results(i) = v.evaluate('metric',obj.metric);

           if obj.verbose
             fprintf(': %f\n',obj.results(i));
           end
           
           if obj.results(i) > maxresult

             maxresult = obj.results(i);
             obj.optimum = i; % optimal configuration
           end

         end
         
         if obj.verbose
           fprintf('optimum configuration ');
           fprintf('%d ',obj.configs(obj.optimum,:));
           fprintf(': %f\n',obj.results(obj.optimum));
         end
         
         % we now know the optimal configuration. Now we retrieve the
         % method of interest with this configuration and retrain
      
         obj.method = obj.validator.procedure.clfmethods{obj.methodidx};

         for j=1:nv
           obj.method.(obj.variables{j}) = obj.configs(obj.optimum,j);
         end
           
         obj = obj.method.train(data,design);
         
       end
       
       function post = test(obj,data)       
           
         post = obj.method.test(data);
         
       end
       
       function p = predict(obj,data)
        % relevant in case the optimizer acts as a predictor
        
        if isa(obj.method,'predictor')
          p = obj.method.predict(data);
        else
          if obj.verbose, fprintf('%s cannot predict\n',class(obj.method)); end
          p = nan(size(data,1),1);
        end
         
       end
       
       function n = nclasses(obj)
        % relevant in case the optimizer acts as a predictor
        
        if isa(obj.method,'predictor')
          n = obj.method.nclasses;
        else
          n = nan;
        end
        
       end

    end
end 
