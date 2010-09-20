classdef ft_mv_ensemble < ft_mv_predictor
%FT_MV_ENSEMBLE takes a number of mvas in parallel and then recombines the
%results using the combfun function
%
% EXAMPLE:
%    m = ft_mv_ensemble('mvas',{ft_mv_naive ft_mv_svm},,'combfun',@(x)(cell2mat(x)))
% creates an ensemble method that combines the results of naive Bayes and
% a support vector machine and puts it in one big array.
%
%   Copyright (c) 2010, Marcel van Gerven


    properties

      mvas              % parallel mva mvas to run
      combfun = @(x)(x) % by default just returns the whole cell array
      parallel = false; % run train function in parallel mode?

    end

    methods
      
       function obj = ft_mv_ensemble(varargin)
         
         obj = obj@ft_mv_predictor(varargin{:});
       
         if isempty(obj.mvas), error('mvas not specified'); end
         
         if ~iscell(obj.mvas), obj.mvas = {obj.mvas}; end
         
         % cast to analysis if necessary
         for c=1:length(obj.mvas)
           if ~isa(obj.mvas{c},'ft_mv_analysis')
             obj.mvas{c} = ft_mv_analysis(obj.mvas{c});
           end
         end
         
       end
       
       function obj = train(obj,X,Y)
     
         if ~iscell(X), X = {X}; end
         if ~iscell(Y), Y = {Y}; end
         
         % we can make one of three choices:
         % - specify multiple datasets and one mva
         % - specify multiple mvas and one dataset 
         % - specify multiple mvas and multiple datasets
         nmva = numel(obj.mvas);
         nsets = length(X);
         if nmva==1 % replicate mvas
           
           if iscell(obj.mvas)
             obj.mvas = repmat(obj.mvas,[1 nsets]);
           else
             obj.mvas = repmat({obj.mvas},[1 nsets]);
           end
           
         else % replicate datasets
           
           if nmva ~= nsets
             X = repmat(X,[1 nmva]);
             Y = repmat(Y,[1 nmva]); 
           end
           
         end
         
         if obj.parallel
         
           if obj.verbose
             fprintf('running %s in parallel mode',class(obj));
           end
           
           obj.mvas = peercellfun(@run_parallel,obj.mvas,repmat({'train'},[1 nmva]),X,Y,'UniformOutput',false);
           
         else
           
           for k=1:length(obj.mvas)
             obj.mvas{k} = obj.mvas{k}.train(X{k},Y{k});
           end
           
         end
         
       end
       
       function Y = test(obj,X)       
           
         Z = cell(1,length(obj.mvas));
         for k=1:length(obj.mvas)
           Z{k} = obj.mvas{k}.test(X);
         end
         
         Y = obj.combfun(Z);
         
       end      
       
       function [m,desc] = model(obj)
         % call the enclosing method
         
         [m,desc] = obj.procedure.model();
           
       end

    end
end 
