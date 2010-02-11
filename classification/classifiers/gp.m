classdef gp < classifier
%GP gaussian process classifier
%
%   Options:
%   'optimize' : if true tries to optimize hyperparameters
%
%   SEE ALSO:
%   gpml-matlab
%
%   Copyright (c) 2008, Marcel van Gerven


    properties
      
      method = 'laplace';
      optimize = true;
      loghyper = [3.0; 0.0];
      
    end
    
    methods
      
      function obj = gp(varargin)
        
        obj = obj@classifier(varargin{:});
        
      end
      
      function p = estimate(obj,X,Y)
        
        if obj.nunique(Y) > 2, error('GP only accepts binary class problems'); end
        
        design = Y;
        
        % make foolproof
        if all(design(:,1) == 1)
          obj.loghyper = -inf;
        elseif all(design(:,1) == 2)
          obj.loghyper = inf;
        else
          
          if ~exist('gpml-matlab','dir')
            error('this code requires an external toolbox: http://www.gaussianprocess.org/gpml/code/matlab/doc/');
          end
          
          targets = design(:,1);
          targets(design == 1) = -1;
          targets(design == 2) = 1;
          
          % training mode
          
          p.data = X;
          p.targets = targets;
          
          if obj.optimize % optimize hyperparameters
            
            if strcmp(obj.method,'laplace')
              
              [p.loghyper p.logmarglik] = minimize(obj.loghyper, 'binaryLaplaceGP', -20, 'covSEiso', 'cumGauss', p.data, p.targets);
              
            else % ep
              
              [p.loghyper p.logmarglik] = minimize(obj.loghyper, 'binaryEPGP', -20, 'covSEiso', p.data, p.targets);
              
            end
          end
         end
      end
      
      function Y = map(obj,X)
        
        if isinf(obj.params.loghyper(1)) 
           if obj.params.loghyper < 0
             Y = [ones(size(X,1),1) zeros(size(X,1),1)];
           else
             Y = [zeros(size(X,1),1) ones(size(X,1),1)];
           end
        else
          
          if strcmp(obj.method,'laplace')
            
            probs = binaryLaplaceGP(obj.params.loghyper, 'covSEiso', 'cumGauss', obj.params.data, obj.params.targets, X);
            
          else % ep
            
            probs = binaryEPGP(obj.params.loghyper, 'covSEiso', obj.params.data, obj.params.targets, X);
            
          end
          
          Y = zeros(size(probs,1),2);
          Y(:,1) = 1 - probs;
          Y(:,2) = probs;
        
        end
        
      end
      
    end
end
