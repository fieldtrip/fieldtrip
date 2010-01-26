classdef circreg < regressor
%CIRCREG circular regression method class
%
%   Circular regression using the Fisher and Lee model.
% 
%   EXAMPLE:
% 
%   effect of mean regressors:
%    c = circreg('mu',0,'beta',[-1 1]','kappa',1);
%    X = randn(1000,2);
%    Y = c.sample(X);
%    scatter(X(:,1),Y));
%    hold on;
%    scatter(X(:,2),Y,'r');
%    xlim([-pi pi]);
%    ylim([-pi pi]);
%    d = c.train(X,Y);
%    fplot(@(x)(d.beta(1)*x),[-pi pi],'b');
%    fplot(@(x)(d.beta(2)*x),[-pi pi],'r');
%    
%   effect of concentration regressors:
%    c = circreg('mu',0,'kappa',1,'verbose',true,'gamma',[1]');
%    X = randn(1000,1);
%    Y = c.sample(X);
%    scatter(X,Y)
%    xlim([-pi pi]);
%    ylim([-pi pi]);
%
%   Copyright (c) 2009, Ali Bahramisharif, Marcel van Gerven
%
%   $Log: circreg.m,v $
%
    properties        
        
        mu            % intercept of von Mises distribution
        beta          % regression coefficients for the mean
        gamma         % regression coefficients for the concentration
        kappa         % concentration parameter of von Mises distribution
                        
        mode =  'mean' % type of circular regression: none, mean, concentration, mixed
              
        lambda = 0;   % regularization parameter for the mean (only defined for gradient descent)
     
    end
    
    methods
        
        function obj = circreg(varargin)
            
            obj = obj@regressor(varargin{:});            
        end        
        
        function obj = train(obj,data,design)
            % determine intercept and regression coefficients                    

            if iscell(data), error('classifier does not take multiple datasets as input'); end
          
            data = data.X;
            design = design.X;
            
            % restrict theta to range [-pi,...,pi]
            rtheta = mod(design(:,1),2*pi);
            rtheta(rtheta > pi) = rtheta(rtheta > pi) - 2*pi;
        
            switch obj.mode
              
              case 'none'
                obj = train_none(obj,data,rtheta);
              case 'mean'
                obj = train_mean(obj,data,rtheta);
              case 'concentration'
                obj = train_concentration(obj,data,rtheta);
              case 'mixed'
                obj = train_mixed(obj,data,rtheta);
              otherwise
                error('unknown mode');
            end
            
        end
        
        function obj = train_none(obj,X,theta)
          % don't regress; just estimate mean and concentration

          % calculate mean
          S = sum(sin(theta));
          C = sum(cos(theta));
          
          obj.mu = angle(1i*S+C);
            
          % calculate concentration
          R = sqrt(C^2 + S^2)/length(theta);          
          if R < 0.53
            obj.kappa = 2*R + R^3 + 5*R^5/6;
          elseif R >= 0.85
            obj.kappa = 1/(R^3 - 4*R^2 + 3*R);
          else
            obj.kappa = -0.4 + 1.39*R + 0.43/(1-R);
          end
          
          if obj.verbose
            disp(obj.loglik(X,theta));
          end
          
        end
        
        function obj = train_mean(obj,X,theta)
          % estimate mean using regressors
          
          if obj.verbose
            fprintf('estimating mean\n');
          end
          
          if isempty(obj.beta) % if not yet initialized
            
            obj.beta = zeros(size(X,2),1);
            
          end
          
          options.Method='lbfgs';
          b = obj.beta;
          obj.beta = minFunc(@(b)regCreg(b(:),theta,X,obj.lambda),b(:),options);
          
          u = X * obj.beta;
          lx = 2*atan(u); % link function
          S = mean(sin(theta - lx));
          C = mean(cos(theta - lx));
          
          % (re)calculate mean and dispersion
          R = sqrt(S^2 + C^2);
          
          % set mean
          obj.mu = angle(1i*S+C);
          
          % set dispersion using A inverse
          if R < 0.53
            obj.kappa = 2*R + R^3 + 5*R^5/6;
          elseif R >= 0.85
            obj.kappa = 1/(R^3 - 4*R^2 + 3*R);
          else
            obj.kappa = -0.4 + 1.39*R + 0.43/(1-R);
          end              
            
        end
        
        function obj = train_concentration(obj,X,theta)
          % estimate concentration using regressors
          
          if obj.verbose
            fprintf('estimating concentration\n');
          end
          
          
          
        end

        function obj = train_mixed(obj,X,theta)
          % estimate mean and concentration using regressors
          
          
        end
        
        function y = loglik(obj,X,theta)
          % log likelihood up to a constant
          
          mus = obj.mu*ones(size(X,1),1);
          if ~isempty(obj.beta)
            mus = mus + 2*atan(X * obj.beta);
          end
          
          kappas = obj.kappa*ones(size(X,1),1);
          if ~isempty(obj.gamma)
            % now obj.kappa plays the role of the offset
            kappas = exp(kappas + X * obj.gamma);
          end
          
         y = -size(X,1) * log(besseli(0,obj.kappa)) + kappas .* sum(cos(theta - mus));
          
        end
        
        function theta = test(obj,data)
          % estimate the angle using intercept and link function
          
          data = data.X;
          
          theta = obj.mu*ones(size(data,1),1);
          if ~isempty(obj.beta)
            theta = theta + 2*atan(data * obj.beta);
          end
          
          theta = dataset(theta);
          
        end
        
        function theta = sample(obj,X)
          % sample from a von Mises distribution (pp. 49 Fisher and Lee)
          
          mus = obj.mu*ones(size(X,1),1);
          if ~isempty(obj.beta)
            mus = mus + 2*atan(X * obj.beta);
          end
          
          kappas = obj.kappa*ones(size(X,1),1);
          if ~isempty(obj.gamma)
            % now obj.kappa plays the role of the offset
            kappas = exp(kappas + X * obj.gamma);
          end
                    
          theta = zeros(size(X,1),1);
          for j=1:length(theta)
            try
              theta(j) = randraw('vonmises',[mus(j) kappas(j)],1);
            catch
              % catch uniform situation
              theta(j) = rand*2*pi;
            end
          end
          
        end
 
    end
end 
