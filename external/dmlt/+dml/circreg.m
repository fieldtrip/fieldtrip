classdef circreg < dml.method
% CIRCREG circular regression method.
%
%   DESCRIPTION
%   Circular regression using the Fisher and Lee model.
%
%   EXAMPLE
%
%   % generate data using fixed mean mu and concentration kappa
%   X = randn(1000,3);
%   c = dml.circreg('mu',0,'kappa',100);
%   Y = c.sample(X);
%   ix = (-pi:0.1:pi)';
%   x = histc(Y,ix); x = x ./ sum(x);
%   polar(ix,1+x);
%
%   % check estimation
%   c = dml.circreg('mode','none','verbose',true);
%   c = c.train([],Y);
%   disp([c.mu c.kappa]);
%
%   % generate data using X-dependent mean mu and concentration kappa
%   X = randn(1000,3);
%   c = dml.circreg('mu',0,'beta',[-1 0 1]','kappa',100);
%   Y = c.sample(X);
%
%   % check estimation
%   c = dml.circreg('mode','mean','verbose',true);
%   c = c.train(X,Y);
%   disp([c.mu c.kappa c.beta']);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)
%   Ali Bahramisharfi (ali@cs.ru.nl)

    properties        
        
        mu            % intercept of von Mises distribution
        kappa         % concentration parameter of von Mises distribution

        beta          % regression coefficients for the mean
        gamma         % regression coefficients for the concentration
        
        method = 1;   % method used to estimate model (1=standard, 2=generalized method of moments, 3=second harmonics)
        
        lambda = 10;  % regularization parameter for the mean (only defined for gradient descent)
        repeat = 1;   % number of repeats for gradient descent (multiple local maxima)
        
        outer = 1;    % maximum number of outer loop iterations in mixed estimation       
        inner = 10;   % maximum number of inner loop iterations in concentration/mixed estimation
        tol = 1;      % smallest update step in concentration/mixed estimation
        
        likelihood    % stored loglikelihood
        
    end
    
    methods
        
        function obj = circreg(varargin)
            
            obj = obj@dml.method(varargin{:});            
        end        
        
        function obj = train(obj,X,Y)
            % determine intercept and regression coefficients                    

            % restrict theta to range [-pi,...,pi]
            rtheta = mod(Y(:,1),2*pi);
            rtheta(rtheta > pi) = rtheta(rtheta > pi) - 2*pi;
        
            if isempty(X), X = nan(size(Y)); end
            
            if isempty(X) || all(isnan(X(:)))
              
              % estimate fixed model
              obj = obj.train_none(X,rtheta);
              
            else
              
              % mean is dependent on X
              obj = obj.train_mean(X,rtheta);
              for j=2:obj.repeat
                o = obj; o.beta = 1e-6*randn(size(X,2),1);
                o = o.train_mean(X,rtheta);
                if o.likelihood > obj.likelihood, obj = o; end
              end
              
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
          
          L = obj.loglik(X,theta);
          if size(X,1)>1 
            obj.likelihood = sum(L);
          end
          if obj.verbose
            fprintf('log likelihood: %.2f\n',obj.likelihood);
          end
          
        end
        
        function obj = train_mean(obj,X,theta)
          % estimate mean using regressors
          
          if obj.verbose, fprintf('estimating mean\n'); end
          
          if isempty(obj.beta) % if not yet initialized
            obj.beta = zeros(size(X,2),1); %1e-6*randn(size(X,2),1);
          end
          
          if ~obj.verbose, options.Display = 'off'; end
          
          b = obj.beta;
          options.MaxIter=1000;
          options.MaxFunEvals=10000;
          b = [0;b];
          if obj.method==1
            obj.beta = minFunc(@(b)regCreg(b(:),theta,X,obj.lambda),b(:),options);
          elseif obj.method==2
            obj.beta = minFunc(@(b)regGMMreg(b(:),theta,X,obj.lambda),b(:),options);
          else
            obj.beta = minFunc(@(b)regCHreg(b(:),theta,X,obj.lambda),b(:),options);
          end
          obj.beta=obj.beta(2:end);
          
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
          
          L = obj.loglik(X,theta);
          if size(X,1)>1 
            obj.likelihood = sum(L);
          end
          if obj.verbose
            fprintf('log likelihood: %.2f\n',obj.likelihood);
          end

        end
        
        function obj = train_concentration(obj,X,theta)
          % estimate concentration using regressors
          % FIXME: Not yet working correctly
          
          if obj.verbose
            fprintf('estimating concentration\n');
          end
          
          obj.gamma = 1e-6*randn(size(X,2),1);
          
          % initialization that escapes local maxima
          nfeatures = size(X,2); nexamples = size(X,1);
          for j=1:nfeatures
            
            [xsrt,xidx] = sort(X(:,j));
            thetasrt = theta(xidx);
            
            m=2; rho = zeros(1,nexamples);
            for i=(1+m):(nexamples-m)
              rho(i) = sqrt(sum(cos(thetasrt((i-m):(i+m))))^2+sum(sin(thetasrt((i-m):(i+m))))^2);
            end
            A1 = besseli(1,rho)./besseli(0,rho);
            A1p = transpose(1 - A1./rho - A1.^2);
            
            %scatter(xsrt((1+m):(nexamples-m)),A1p((1+m):(nexamples-m)));
            model = regress(A1p,[ones(size(X,1),1) X(:,j)]);
            obj.kappa = model(1);
            obj.gamma(j) = model(2);
            
          end
          
          X1 = [ones(size(X,1),1) X];
          
          for idx2=0:obj.inner
            
            oldgamma = obj.gamma;
            
            bk = exp(obj.kappa + X * obj.gamma);
            
            % compute mean
            S = bk' * sin(theta);
            C = bk' * cos(theta);
            R = sqrt(S^2 + C^2);
            obj.mu = angle(1i*S+C);
            
            % update equation for regression coefficients
            % G2 = W
            
            % compute W
            A1 = besseli(1,bk)./besseli(0,bk);
            A1(isinf(A1)) = 1; % bug fix: check validity
            
            A1p = 1 - A1./bk - A1.^2;
            W = diag(bk.^2 .*  A1p);
            
            % compute y
            by = (cos(theta - obj.mu) - A1) ./ (bk .*  A1p);
            
            WG2 = X1' * W;
            
            delta = inv(WG2 * X1) * (WG2 * by);
            if ~any(isnan(delta))
              obj.kappa = obj.kappa + delta(1);
              obj.gamma = delta(2:end) + obj.gamma;
            end
            
            % jump out in case of very small updates
            if norm(obj.gamma - oldgamma,2) < obj.tol, break; end
            
            L = obj.loglik(X,theta);
            if size(X,1)>1
              obj.likelihood = [obj.likelihood sum(L)];
            end
            if obj.verbose
              fprintf('log likelihood: %.2f\n',obj.likelihood(end));
            end

          end
          
        end

        function obj = train_mixed(obj,X,theta)
          % estimate mean and concentration using regressors
          % FIXME: not yet working correctly
          
          % find starting values
          obj = train_mean(obj,X,theta);
          
          obj.gamma = 1e-6*randn(size(X,2),1);
          
          % initialization that escapes local maxima
          nfeatures = size(X,2); nexamples = size(X,1);
          for j=1:nfeatures
            
            [xsrt,xidx] = sort(X(:,j));
            thetasrt = theta(xidx);
            
            m=2; rho = zeros(1,nexamples);
            for i=(1+m):(nexamples-m)
              rho(i) = sqrt(sum(cos(thetasrt((i-m):(i+m))))^2+sum(sin(thetasrt((i-m):(i+m))))^2);
            end
            A1 = besseli(1,rho)./besseli(0,rho);
            A1p = transpose(1 - A1./rho - A1.^2);
            
            model = regress(A1p,[ones(size(X,1),1) X(:,j)]);
            obj.kappa = model(1);
            obj.gamma(j) = model(2);
          end
          
          if obj.verbose
            fprintf('estimating mixed model\n');
          end
          
          X1 = [ones(size(X,1),1) X];
          
          for idx1=0:obj.outer
            
            oldbeta1 = obj.beta;
            oldgamma1 = obj.gamma;
            
            if ~isempty(obj.gamma)
              u = X * obj.gamma;
            else
              u = 0;
            end
            bk = exp(obj.kappa + u);
            A1 = besseli(1,bk)./besseli(0,bk);
            A1(isinf(A1)) = 1; % bug fix: check validity
            
            % update betas
            for idx2=1:obj.inner
              
              oldbeta2 = obj.beta;
              
              v = X * obj.beta;
              
              bu = sin(theta - obj.mu - 2*atan(v));
              dlx = 2./(1+v.^2); % derivative of link function
              G = diag(dlx);
              by = bu./(A1.*dlx);
              
              K = diag(bk .* A1);
              WG2 = X' * G * K * G;
              
              % update equation for regression coefficients
              obj.beta = inv(WG2 * X) * (WG2 * by) + obj.beta;
              
              % jump out in case of very small updates
              if norm(obj.beta - oldbeta2,2) < obj.tol, break; end
              
            end
            
            % (re)calculate mean
            lx = 2*atan(X * obj.beta); % link function
            S = mean(sin(theta - lx));
            C = mean(cos(theta - lx));
            R = sqrt(S^2 + C^2);
            obj.mu = angle(1i*S+C);
            
            if ~isempty(obj.beta)
              v = 2*atan(X * obj.beta);
            else
              v = 0;
            end
            
            % update kappas
            for idx2=0:obj.inner
              
              oldgamma2 = obj.gamma;
              
              if ~isempty(obj.gamma)
                u = X * obj.gamma;
              else
                u = 0;
              end
              bk = exp(obj.kappa + u);
              
              % compute mean
              S = bk' * sin(theta);
              C = bk' * cos(theta);
              R = sqrt(S^2 + C^2);
              %       obj.mu = asin(S/R);
              
              % update equation for regression coefficients
              % G2 = W
              
              % compute W
              A1 = besseli(1,bk)./besseli(0,bk);
              
              % bug fix: check validity
              A1(isinf(A1)) = 1;
              
              A1p = 1 - A1./bk - A1.^2;
              W = diag(bk.^2 .*  A1p);
              
              % compute y
              by = (cos(theta - obj.mu - v) - A1) ./ (bk .*  A1p);
              
              WG2 = X1' * W;
              
              delta = inv(WG2 * X1) * (WG2 * by);
              if ~any(isnan(delta))
                obj.kappa = obj.kappa + delta(1);
                obj.gamma = obj.gamma + delta(2:end);
              end
              
              % jump out in case of very small updates
              if norm(obj.gamma - oldgamma2,2) < obj.tol, break; end
              
            end
            
            % jump out in case of very small updates
            if norm(obj.beta - oldbeta1,2) < obj.tol && ...
                norm(obj.gamma - oldgamma1,2) < obj.tol
              break;
            end
            
            L = obj.loglik(X,theta);
            if size(X,1)>1
              obj.likelihood = [obj.likelihood sum(L)];
            end
            if obj.verbose
              fprintf('log likelihood: %.2f\n',obj.likelihood(end));
            end
            
          end
          
        end
        
        function theta = test(obj,X)
          % estimate the angle using intercept and link function
          
          theta = obj.mu*ones(size(X,1),1);
          if ~isempty(obj.beta)
            theta = theta + 2*atan(X * obj.beta);
          end
           
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
        
        function rtheta = sample(obj,X)
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
          
          % map the thing back to [-pi pi]
          rtheta = mod(theta(:,1),2*pi);
          rtheta(rtheta > pi) = rtheta(rtheta > pi) - 2*pi;
          
        end
 
    end
end 