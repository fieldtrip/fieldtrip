classdef ft_mv_analysis
%FT_MV_ANALYSIS multivariate analysis class
%   
%   A multivariate analysis is just a sequence of multivariate methods 
%   {method1 method2 method3 ...} that are called in this order and where 
%   the output of the previous method acts
%   as input to the next method.
%
%   Copyright (c) 2008, Marcel van Gerven

    properties
      
      mvmethods; % the methods that specify the classification procedure
        
    end
    
    methods
      
       function obj = ft_mv_analysis(mvmethods)
       % constructor expects mva methods

       if ~nargin, error('methods not specified'); end
        
        if isa(mvmethods,'ft_mv_analysis')
          obj = mvmethods;
          return;
        end
        
        % cast to cell array if only one method is specified
        if ~iscell(mvmethods), mvmethods = { mvmethods }; end
        
        obj.mvmethods = mvmethods;
       
       end     
       
       function obj = train(obj,X,Y)
         % train just calls the methods' train functions in order to
         % produce an output
         
         if nargin<3, Y = []; end % preprocessors don't expect a design
      
         % transform input to matrix
         if iscell(X)
           for c=1:length(X)
              if ndims(X{c})~=2
                sz = size(X{c});
                X{c} = reshape(X{c},[sz(1) prod(sz(2:end))]);
              end
           end
         else
           if ndims(X)~=2
             sz = size(X);
             X = reshape(X,[sz(1) prod(sz(2:end))]);
           end
         end
         
         % transform output to matrix
         if iscell(Y)
           for c=1:length(Y)
              if ndims(Y{c})~=2
                sz = size(Y{c});
                Y{c} = reshape(Y{c},[sz(1) prod(sz(2:end))]);
              end
           end
         else
           if ndims(Y)~=2
             sz = size(Y);
             Y = reshape(Y,[sz(1) prod(sz(2:end))]);
           end
         end
         
         for c=1:length(obj.mvmethods)

           obj.mvmethods{c} = obj.mvmethods{c}.train(X,Y);
           
           if c<length(obj.mvmethods)
             X = obj.mvmethods{c}.test(X);
           end
             
         end
         
       end
       
       function Y = test(obj,X)
         % test just calls the methods' test functions in order to produce a posterior
         
         if isempty(X)
           Y = [];
           return;
         end
         
         Y = X;

          % transform output to matrix
         if iscell(Y)
           for c=1:length(Y)
              if ndims(Y{c})~=2
                sz = size(Y{c});
                Y{c} = reshape(Y{c},[sz(1) prod(sz(2:end))]);
              end
           end
         else
           if ndims(Y)~=2
             sz = size(Y);
             Y = reshape(Y,[sz(1) prod(sz(2:end))]);
           end
         end
         
         for c=1:length(obj.mvmethods)
           Y = obj.mvmethods{c}.test(Y);
         end
         
       end
       
       function s = name(obj)
         % returns multivariate analysis as a string
         
         s = subname(obj.mvmethods);
         
         function ss = subname(m)
            
           if iscell(m)
             ss = ['{ ' cell2mat(cellfun(@(x)([subname(x) ' ']),m,'UniformOutput',false)) '}'];
           else
             ss = class(m);
           end
           
         end
         
       end
       
       function [m,desc] = model(obj)
         % just get the final model implied by this mva
       
           [m,desc] = obj.mvmethods{end}.model();
         
       end
       
    end
    
end
