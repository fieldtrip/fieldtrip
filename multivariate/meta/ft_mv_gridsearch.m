classdef ft_mv_gridsearch < ft_mv_meta
%FT_MV_GRIDSEARCH can be used to optimize certain parameters of a multivariate analysis
%
%   EXAMPLES:
%  
%   ft_mv_gridsearch('verbose',true,'mva',ft_mv_svm,'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'vars','C','vals',logspace(-3,3,7))
%
%     will act as an svm but will optimize the variable
%     C in the range logspace(-3,3,7) using the specified validator.
%
%   ft_mv_gridsearch('verbose',true,'mva',{ft_mv_filterer ft_mv_svm},'validator',ft_mv_crossvalidator('nfolds',0.8,'metric','accuracy'),'mvidx',[1 2],'vars',...
%         {'maxfeatures' 'C'},'vals',{1:10 [1 10]})
%
%     will optimize the combined assignments of both the filterer and the svm
%
%   NOTE:
%     mvmethods are supposed to handle optimization efficiently by starting
%     at previous solutions
%
%   Copyright (c) 2009, Marcel van Gerven


    properties

      parallel = false % parallel processing
      
      validator % the used validator
      
      mvidx     % the index of the mvmethod in the mva that is to be optimized
      vars      % the variables to optimize
      vals      % the values used
      
      configs   % the configurations (cartesian product of vals)
      results   % result per configuration
      
      optimum   % optimal configuration
      
    end

    methods
      
       function obj = ft_mv_gridsearch(varargin)
         
          obj = obj@ft_mv_meta(varargin{:});
         
         if isempty(obj.validator)
           obj.validator = ft_mv_crossvalidator('metric','accuracy');
         end
         
         assert(~isempty(obj.mva));
         assert(~isempty(obj.vars));
         assert(~isempty(obj.vals));
         
         obj.mva = ft_mv_analysis(obj.mva);
          
         if isempty(obj.mvidx)
           if length(obj.mva.mvmethods) == 1
             obj.mvidx = ones(1,length(obj.vars));
           else
             error('mvidx not specified');
           end
         end

         if ~iscell(obj.vars)
           obj.vars = { obj.vars };
           obj.vals = { obj.vals };
         end
         
         obj.validator.mva = obj.mva;
         
       end
       
       function obj = train(obj,X,Y)
         
         if obj.verbose
           fprintf('optimizing variables');
           for c=1:length(obj.vars)
             fprintf('%s.%s',class(obj.mva.mvmethods{obj.mvidx(c)}),obj.vars{c});
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
         
         if obj.parallel
           train_parallel();
         else
           train_sequential();
         end
                  
         if obj.verbose
           fprintf('optimum configuration ');
           for j=1:size(obj.configs,2)
               fprintf(' %s.%s=%f',class(obj.mva.mvmethods{obj.mvidx(j)}),obj.vars{j},obj.configs(obj.optimum,j));
           end
           fprintf(' : %f\n',obj.results(obj.optimum));
         end
         
         % we now know the optimal configuration. Now we retrieve the
         % method of interest with this configuration and retrain with all
         % data
         
         if isempty(obj.optimum)
           error('no optimum configuration found; check performance measure');
         end
      
         for j=1:nv
           obj.mva.mvmethods{obj.mvidx(j)}.(obj.vars{j}) = obj.configs(obj.optimum,j);
         end
         obj.mva = obj.mva.train(X,Y);
         
         function train_sequential()
           
           vld = obj.validator; % the mva used by the validator
             
           for i=1:size(obj.configs,1)
             
             % set parameters
             if iscell(vld.mva)
               % happens when the vld has already been validated
               for c=1:numel(vld.mva)
                 for j=1:nv
                   vld.mva{c}.mvmethods{obj.mvidx(j)}.(obj.vars{j}) = obj.configs(i,j);
                 end
               end
             else
               for j=1:nv
                 vld.mva.mvmethods{obj.mvidx(j)}.(obj.vars{j}) = obj.configs(i,j);
               end
             end
             
             if obj.verbose
               fprintf('evaluating configuration ');
               for c=1:length(obj.vars)
                 fprintf(' %s.%s=%f',class(obj.mva.mvmethods{obj.mvidx(c)}),obj.vars{c},obj.configs(i,c));
               end
               fprintf('\n');
             end
             
             vld = vld.train(X,Y);
             
             obj.results(i) = vld.performance;
             
             if obj.verbose
               fprintf('%s : %f\n',vld.metric,obj.results(i));
             end
             
             if obj.results(i) > maxresult
               
               maxresult = obj.results(i);
               obj.optimum = i; % optimal configuration
             end
             
           end
           
         end
         
         function train_parallel()
           
           % create all validators
           nconf = size(obj.configs,1);
           vld = cell(nconf,1);
           for i=1:nconf
             
             vld{i} = obj.validator; % the mva used by the validator
             
             % set parameters
             if iscell(vld{i}.mva)
               % happens when the vld has already been validated
               for c=1:numel(vld.mva)
                 for j=1:nv
                   vld{i}.mva{c}.mvmethods{obj.mvidx(j)}.(obj.vars{j}) = obj.configs(i,j);
                 end
               end
             else
               for j=1:nv
                 vld{i}.mva.mvmethods{obj.mvidx(j)}.(obj.vars{j}) = obj.configs(i,j);
               end
             end
             
           end
             
           vld = peercellfun(@run_parallel,vld,repmat({'train'},[1 nconf]),repmat({X},[1 nconf]),repmat({Y},[1 nconf]),'UniformOutput',false);
         
           for i=1:nconf
             obj.results(i) = vld{i}.performance;
             
             if obj.results(i) > maxresult
               maxresult = obj.results(i);
               obj.optimum = i; % optimal configuration
             end
             
           end
           
           
         end
         
       end
       
       function Y = test(obj,X)       
           
         Y = obj.mva.test(X);
         
       end      
       
       function [m,desc] = model(obj)
         % call the enclosing method
         
         [m,desc] = obj.mva.model();
           
       end

    end
end 
