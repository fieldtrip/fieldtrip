classdef linreg < regressor
%LINREG linear regression method class
%
% Note:
% ridge regression allows multiple outputs
% elastic net allows L2 to be a structure matrix
%
% Copyright (c) 2008, Marcel van Gerven

  properties
    
    L1 = 0; % L1 regularization parameter
    L2 = 0; % L2 regularization parameter
  
  end
    
    methods
        
        function obj = linreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function p = estimate(obj,X,Y)
          
          % handle zero vector
          if ~any(Y)
            p.model = zeros(size(X,2)+1,1);
            return;
          end
          
          if obj.L1 == 0 && all(obj.L2(:) == 0)
            % unregularized 
            
            if obj.verbose
              fprintf('computing unregularized linear regression\n');
            end
            
            p.model = regress(Y,[X ones(size(X,1),1)]); % X \ Y
          
          elseif obj.L1 == 0 && all(obj.L2(:) ~= 0)
            % ridge regression
            
            if obj.verbose
              fprintf('computing ridge regression with L2 parameter %g\n',obj.L2);
            end

            % penalize the absolute value of each element by the same amount
            % don't penalize bias term 
            lambdas = [obj.L2*ones(size(X,2),1); 0];

            X = [X ones(size(X,1),1)];
            R = chol(X'*X + diag(lambdas));
            
            p.model = R\(R'\(X'*Y));
            
          elseif obj.L1 ~= 0 && all(obj.L2(:) == 0)
            
            if obj.verbose
              fprintf('computing lasso with L1 parameter %g\n',obj.L1);
            end
            
            % penalize the absolute value of each element by the same amount
            % don't penalize bias term
            lambdas = [obj.L1*ones(size(X,2),1); 0];
            
            X = [X ones(size(X,1),1)];
            
            R = chol(X'*X + diag(lambdas));
            
            funObj = @(w)GaussianLoss(w,X,Y); % Loss function that L1 regularization is applied to
            
            if ~isempty(obj.params) && isfield(obj.params,'model') && ~isempty(obj.params.model)
              % warm start
              w_init = obj.params.model;
            else
              w_init = R\(R'\(X'*Y)); % Initial value for iterative optimizer
            end
            
            p.model = L1GeneralProjection(funObj,w_init,lambdas);
            
          else % elastic net
            
            if ~isempty(obj.params) && isfield(obj.params,'model') && ~isempty(obj.params.model)
              % warm start
              w_init = obj.params.model;
              [beta,beta0] = elastic(X',Y',obj.L1,obj.L2,[],w_init(1:(end-1)),w_init(end));
            else
              [beta,beta0] = elastic(X',Y',obj.L1,obj.L2);
            end
            
            p.model = [beta; beta0];

          end
          
          if ~any(p.model)
            warning('zero vector as model; too strong regularization?');
          end
          
        end
        
        function Y = map(obj,X)
          Y = [X ones(size(X,1),1)] * obj.params.model;
        end
        
        function [m,desc] = getmodel(obj)
          % return the parameters in some shape
          
          m = {obj.params.model(1:(end-1),:)'}; % ignore bias term
          desc = {'unknown'};
          
        end

    end
end 
