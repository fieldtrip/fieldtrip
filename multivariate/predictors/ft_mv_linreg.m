classdef ft_mv_linreg < ft_mv_predictor
%FT_MV_LINREG linear regression method class
%
% Note:
% ridge regression allows multiple outputs
% elastic net allows L2 to be a structure matrix
%
% Copyright (c) 2008, Marcel van Gerven

  properties
    
    L1 = 0; % L1 regularization parameter
    L2 = 0; % L2 regularization parameter
    
    weights; % linear regression coefficients
  
  end
    
    methods
        
        function obj = ft_mv_linreg(varargin)
            
            obj = obj@ft_mv_predictor(varargin{:});            
        end        
        
        function obj = train(obj,X,Y)
          
          % multiple datasets
          if iscell(X) || iscell(Y)
            obj = ft_mv_ndata('mvmethod',obj);
            obj = obj.train(X,Y);
            return;
          end
          
          % missing data
          if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
          
          % multiple outputs
          if size(Y,2) > 1
            obj = ft_mv_noutput('mvmethod',obj);
            obj = obj.train(X,Y);
            return;
          end
          
          % handle zero vector
          if ~any(Y)
            obj.weights = zeros(size(X,2)+1,1);
            return;
          end
          
          if obj.L1 == 0 && all(obj.L2(:) == 0)
            % unregularized 
            
            if obj.verbose
              fprintf('computing unregularized linear regression\n');
            end
            
            obj.weights = regress(Y,[X ones(size(X,1),1)]); % X \ Y
          
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
            
            obj.weights = R\(R'\(X'*Y));
            
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
            
            if ~isempty(obj.params) && isfield(obj.params,'weights') && ~isempty(obj.params.weights)
              % warm start
              w_init = obj.params.weights;
            else
              w_init = R\(R'\(X'*Y)); % Initial value for iterative optimizer
            end
            
            opt.verbose = 0;
            obj.weights = L1GeneralProjection(funObj,w_init,lambdas,opt);
            
          else % elastic net
            
            if ~isempty(obj.params) && isfield(obj.params,'weights') && ~isempty(obj.params.weights)
              % warm start
              w_init = obj.params.weights;
              [beta,beta0] = elastic(X',Y',obj.L1,obj.L2,[],w_init(1:(end-1)),w_init(end));
            else
              [beta,beta0] = elastic(X',Y',obj.L1,obj.L2);
            end
            
            obj.weights = [beta; beta0];

          end
          
          if ~any(obj.weights)
            warning('zero vector as weights; too strong regularization?');
          end
          
        end
        
        function Y = test(obj,X)
          Y = [X ones(size(X,1),1)] * obj.weights;
        end
        
        function [m,desc] = model(obj)
          % return the parameters in some shape
          
          m = {obj.weights(1:(end-1),:)'}; % ignore bias term
          desc = {'linear regression coefficients'};
          
        end

    end
end 
