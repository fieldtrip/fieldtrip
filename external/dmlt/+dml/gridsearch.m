classdef gridsearch < dml.method
% GRIDSEARCH grid search method.
%
%   DESCRIPTION
%   This method can be used to optimize certain parameters of a
%   multivariate analysis as specified in a crossvalidator object. The
%   method will return an optimized multivariate method and hence can be
%   used in a more complex multivariate analysis pipeline.
%   Note: methods are supposed to handle optimization efficiently by warm starting
%   at previous solutions. This also requires the 'restart' parameter of
%   this method to be set to false. This is done to prevent unwanted use of
%   previously estimated parameters by an object whenever 'restart' is true.
%   
%
%   EXAMPLE:
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]';
%   v = dml.enet.lambdapath(X,Y,'logistic',5,1e-2);
%   m = dml.gridsearch('validator',dml.crossvalidator('type','split','stat','accuracy','mva',dml.enet('type','logistic','restart',false)),'vars','L1','vals',v,'verbose',true);
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)


    properties

      validator % the used cross-validator
      
      idx       % the index of the method in the mva that is to be optimized
      vars      % the variables to optimize
      vals      % the values used
      
      configs   % the configurations (cartesian product of vals)
      outcome   % outcome per configuration
      
      optimum   % optimal configuration
      
      mva       % optimized mva method
      
      models    % all models learned for each of the parameter settings
      
      retrain = true   % retrain the model using the optimal setting on all data
      
    end

    methods
      
       function obj = gridsearch(varargin)
         
	     obj = obj@dml.method(varargin{:});

         assert(~isempty(obj.validator));
         assert(~isempty(obj.vars));
         assert(~isempty(obj.vals));
         
         obj.mva = obj.validator.mva;
         if ~isa(obj.mva,'dml.analysis')
          obj.mva = dml.analysis(obj.mva);
         end
         
         if ~iscell(obj.vars)
           obj.vars = { obj.vars };
           obj.vals = { obj.vals };
         end

         if isempty(obj.idx)
           obj.idx = length(obj.mva.method)*ones(1,length(obj.vars));
         end
         
         if ~obj.retrain
           assert(~obj.validator.compact); % we should remember trained objects
         end
         
       end
       
       function obj = train(obj,X,Y)
                  
         if obj.verbose
           fprintf('optimizing variable(s)');
           for c=1:length(obj.vars)
             fprintf(' %s.%s',class(obj.mva.method{obj.idx(c)}),obj.vars{c});
           end
           fprintf('\n');
         end
         
         % the optimizer iterates over all values for the specified
         % variables and returns the best result
         
         nv = length(obj.vars);
         
         maxresult = -Inf;
         
         % iterate over configurations
         if length(obj.vals)==1
           % allows ordering from high to low
           % can be faster in some cases with warm restarts
           obj.configs = obj.vals{1}(:);
         else
           obj.configs = cartprod(obj.vals{:});
         end
         
         maxvld = [];
         train_sequential();
         
         if isempty(obj.optimum)
           error('no optimum configuration found; check performance measure');
         end
         
         if obj.retrain

           if obj.verbose
             fprintf('retraining optimum configuration');
             for j=1:size(obj.configs,2)
               fprintf(' %s.%s=%f',class(obj.mva.method{obj.idx(j)}),obj.vars{j},obj.configs(obj.optimum,j));
             end
             fprintf(' using all data\n');
           end
           
           % we now know the optimal configuration. Now we retrieve the
           % method of interest with this configuration and retrain with all
           % data
           
           for j=1:nv
             obj.mva.method{obj.idx(j)}.(obj.vars{j}) = obj.configs(obj.optimum,j);
           end
           obj.mva = obj.mva.train(X,Y);
           
         else
           
          obj.mva = maxvld.mva{1};
          maxvld = [];
           
         end
           
         function train_sequential()
           
           vld = obj.validator;
           
           obj.outcome = zeros(1,size(obj.configs,1));
           obj.models = cell(1,size(obj.configs,1));
           for i=1:size(obj.configs,1)
             
             % set parameters
             if iscell(vld.mva)
               % happens when the vld has already been validated
               for c=1:numel(vld.mva)
                 for j=1:nv
                   vld.mva{c}.method{obj.idx(j)}.(obj.vars{j}) = obj.configs(i,j);
                 end
               end
             else
               for j=1:nv
                 vld.mva.method{obj.idx(j)}.(obj.vars{j}) = obj.configs(i,j);
               end
             end
             
             if obj.verbose
               fprintf('evaluating configuration %d of %d:',i,size(obj.configs,1));
               for c=1:length(obj.vars)
                 fprintf(' %s.%s=%f',class(obj.mva.method{obj.idx(c)}),obj.vars{c},obj.configs(i,c));
               end
               fprintf('\n');
             end
             
             vld = vld.train(X,Y);
             
             obj.outcome(i) = mean(vld.statistic);
             
             obj.models{i} = vld.model;
             
             if obj.verbose
               fprintf('%s : %f\n',vld.statistic,obj.outcome(i));
             end
             
             if obj.outcome(i) > maxresult
               
               maxresult = obj.outcome(i);
               if ~obj.retrain
                 maxvld = vld;
               end
               obj.optimum = i; % optimal configuration
             end
             
           end
           
         end
         
       end
         
       function Y = test(obj,X)
         
         Y = obj.mva.test(X);
         
       end
       
       function m = model(obj)
         % call the enclosed method
         
         m = obj.mva.model();
         
       end
       
    end
    
end
