classdef ft_mv_ensemble < ft_mv_meta
%FT_MV_ENSEMBLE takes a number of mva in parallel and then recombines the
%results using the combfun function
%
% EXAMPLE:
%    m = ft_mv_ensemble('mva',{ft_mv_naive ft_mv_svm},'combfun',@(x)(cell2mat(x)))
% creates an ensemble method that combines the results of naive Bayes and
% a support vector machine and puts it in one big array.
%
%   Copyright (c) 2010, Marcel van Gerven


    properties

      combfun = @(x)(x); % by default just returns the whole cell array
      parallel = false;  % run train function in parallel mode?

    end

    methods
      
       function obj = ft_mv_ensemble(varargin)
         
         obj = obj@ft_mv_meta(varargin{:});
       
         if isempty(obj.mva), error('mva not specified'); end
         
         if ~iscell(obj.mva), obj.mva = {obj.mva}; end
         
         if iscell(obj.mva)
           for c=1:length(obj.mva)
             if ~isa(obj.mva{c},'ft_mv_analysis')
               obj.mva{c} = ft_mv_analysis(obj.mva{c});
             end
           end
         else
           if ~isa(obj.mva,'ft_mv_analysis')
             obj.mva = ft_mv_analysis(obj.mva);
           end
         end
         
       end
       
       function obj = train(obj,X,Y)
     
         if ~iscell(X), X = {X}; end
         if ~iscell(Y), Y = {Y}; end
         
         % we can make one of three choices:
         % - specify multiple datasets and one mva
         % - specify multiple mva and one dataset 
         % - specify multiple mva and multiple datasets
         nmva = numel(obj.mva);
         nsets = length(X);
         if nmva==1 % replicate mva
           
           if iscell(obj.mva)
             obj.mva = repmat(obj.mva,[1 nsets]);
           else
             obj.mva = repmat({obj.mva},[1 nsets]);
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
           
           obj.mva = peercellfun(@run_parallel,obj.mva,repmat({'train'},[1 nmva]),X,Y,'UniformOutput',false);
           
         else
           
           for k=1:length(obj.mva)
             obj.mva{k} = obj.mva{k}.train(X{k},Y{k});
           end
           
         end
         
       end
       
       function Y = test(obj,X)       
           
         Z = cell(1,length(obj.mva));
         for k=1:length(obj.mva)
           Z{k} = obj.mva{k}.test(X);
         end
         
         Y = obj.combfun(Z);
         
       end      
       
       function [m,desc] = model(obj)
         % call the enclosing method
         
         [m,desc] = obj.mva.model();
           
       end

    end
end 
