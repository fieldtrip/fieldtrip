classdef mvmethod
% MVMETHOD base class for multivariate methods
%   
%   This base class contains common properties
%   which may be called by all child methods
%
%   mainly deals with data handling
%
%   Copyright (c) 2009, Marcel van Gerven


  properties
  
    verbose = false;
    
    params; % the used parameters for the mapping/unmapping
    
    % here we store some useful properties of training data and design
    
    indims    % data dimensions
    outdims   % design dimensions
    
  end
  
  methods
    
    function [m,desc] = getmodel(obj)
      % default behaviour when we ask for a model (override in subclass)
      
      if obj.verbose
        fprintf('don`t know how to return model for object of type %s; returning empty model and description\n',class(obj));
      end
      
      m = {};
      desc = {};
      
    end
    
    function obj = train(obj,data,design)
      
      if iscell(data) && ~obj.istransfer()
        % compute result for each individual dataset
        
        params = cell(1,length(data));
        
        for c=1:length(data)
          obj = obj.train(data{c},design{c});
          params{c} = obj.params;
        end
        
        obj.params = params;
        
      else
        
        bindim = isempty(obj.indims);
        boutdim = isempty(obj.outdims);
        
        % data and design are collapsed to matrices
        if iscell(data) && obj.istransfer()
          
          if length(unique(cellfun(@(x)(size(x,2)),data))) > 1
            % check if datasets have the same number of features
            error('datasets must have the same number of features for transfer learning');
          end
          
          if bindim
            obj.indims = cell(1,length(data));
          end
          for c=1:length(data)
            if bindim
              obj.indims{c} = size(data{c});
              obj.indims{c} = obj.indims{c}(2:end);
            end
            data{c} = data{c}(1:size(data{c},1),:);
          end
          
          if boutdim
            obj.outdims = cell(1,length(design));
          end
          for c=1:length(design)
            if boutdim
              obj.outdims{c} = size(design{c});
              obj.outdims{c} = obj.outdims{c}(2:end);
            end
            design{c} = design{c}(1:size(design{c},1),:);
          end
          
        else
          
          if bindim
            obj.indims = size(data);
            obj.indims = obj.indims(2:end);
          end
          
          if ndims(data)~=2
            data = data(1:size(data,1),:);
          end
          
          if boutdim
            obj.outdims = size(design);
            obj.outdims = obj.outdims(2:end);
          end
          
          if ndims(design)~=2
            design = design(1:size(design,1),:);
          end
          
        end
        
        obj.params = obj.estimate(data,design);
        
      end
    end
    
    function data = test(obj,data)
      
      if iscell(data) && ~obj.istransfer()
        
        p = obj.params;        
        for c=1:length(data)
          
          obj.params = p{c};
          data{c} = obj.test(data{c});
        end
        
      else
        
        % data is collapsed to a matrix
        if iscell(data)
          for c=1:length(data)
            if ndims(data{c})~=2
              data{c} = data{c}(1:size(data{c},1),:);
            end
          end
        else
          if ndims(data)~=2
            data = data(1:size(data,1),:);
          end
        end
        
        data = obj.map(data);
        
        % try to map result back to original dimensions
        
        if iscell(data)
          for c=1:length(data)
            if numel(data{c}) == prod([size(data{c},1) obj.indims{c}])
              data{c} = reshape(data{c},[size(data{c},1) obj.indims{c}]);
            end
          end
        else
          
          % this fix is needed when a combiner maps multiple datasets
          % into one dataset
          if iscell(obj.indims)
            obj.indims = obj.indims{1};
          end
          
          if numel(data) == prod([size(data,1) obj.indims])
            data = reshape(data,[size(data,1) obj.indims]);
          end
        end
        
      end
      
    end
    
    function data = untest(obj,data)
      % invert the mapping
      
      if iscell(data) && ~obj.istransfer()
        
        p = obj.params;
        
        for c=1:length(data)
          
          obj.params = p{c};
          data{c} = obj.untest(data{c});
        end
        
      else
        
        % data is collapsed to a matrix
        if iscell(data)
          for c=1:length(data)
            if ndims(data{c})~=2
              data{c} = data{c}(1:size(data{c},1),:);
            end
          end
        else
          if ndims(data)~=2
            data = data(1:size(data,1),:);
          end
        end
        
        data = obj.unmap(data);
        
        % try to map result back to original dimensions
        
        if iscell(data)
          for c=1:length(data)
            if numel(data{c}) == prod([size(data{c},1) obj.indims{c}])
              data{c} = reshape(data{c},[size(data{c},1) obj.indims{c}]);
            end
          end
        else
          if numel(data) == prod([size(data,1) obj.indims])
            data = reshape(data,[size(data,1) obj.indims]);
          end
        end
        
      end
    end    
    
    function p = estimate(obj,X,Y)
      % parameter estimation
      p = [];
      
    end
    
    function Y = map(obj,X)
      % default identity mapping
      
      Y = X;
      
    end
    
    function X = unmap(obj,Y)
      % default identity mapping
      
      X = Y;
      
    end
    
    function b = istransfer(obj)
      % return whether or not this method is a transfer learner
      % must be overloaded by e.g., one_against_one
      
      if isa(obj,'transfer_learner')
        b = true;
      else
        b = false;
      end
    end
    
  end
  
  
  
  methods(Static=true)
    % some helper functions operating on datasets
    
    function Y = labeled(X)
      % return indices of labeled (non-nan) datapoints
      Y = find(any(~isnan(X(1:size(X,1),:)),2));
    end
    
    function Y = unlabeled(X)
      % return indices of unlabeled (nan) datapoints
      Y = find(any(isnan(X(1:size(X,1),:)),2));
    end
    
    function Y = unique(X)
      % return the unique trials
      Y = unique(X(1:size(X,1),:),'rows');
    end
    
    function n = nunique(X)
      % return the number of unique trials
      [tmp,tmp,idx] = unique(X(1:size(X,1),:),'rows');
      n = max(idx);
    end
    
  end
  
end