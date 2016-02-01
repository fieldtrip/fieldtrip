function waic = gp_waic(gp, x, y, varargin)
%GP_WAIC The widely applicable information criterion (WAIC) for GP model
% 
%  Description
%    WAIC = GP_WAIC(GP, X, Y) evaluates WAIC defined by
%    Watanabe(2010) given a Gaussian process model GP, training
%    inputs X and training outputs Y. Instead of Bayes loss we
%    compute the Bayes utility which is just the negative of loss
%    used by Watanabe.
% 
%    WAIC is evaluated as follows when using the variance form
%        
%      WAIC(n) = BUt(n) - V/n
%        
%    where BUt(n) is Bayesian training utility,  V is functional variance
%    and n is the number of training inputs.
%
%      BUt = mean(log(p(yt | xt, x, y)))
%      V = sum(E[log(p(y|th))^2] - E[log(p(y|th))]^2)
%
%    When using the Gibbs training loss, WAIC is evaluated as follows
%
%          WAIC(n) = BUt(n) - 2*(BUt(n) - GUt(n))
%
%    where BUt(n) is as above and GUt is Gibbs training utility
%
%          GUt(n) = E_th[mean(log(p(y|th)))].
%     
%    GP can be a Gaussian process structure, a record structure
%    from GP_MC or an array of GPs from GP_IA.
%
%   OPTIONS is optional parameter-value pair
%      method - Method to evaluate waic, 'V' = Variance method, 'G' = Gibbs
%               training utility method (default = 'V')
%      form -   Return form, 'mean' returns the mean value and 'all'
%               returns the values for all data points (default = 'mean')
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%
%   See also
%     GP_DIC, DEMO_MODELASSESMENT1, DEMO_MODELASSESMENT2
%
%   References
%     
%     Watanabe(2010). Equations of states in singular statistical
%     estimation. Neural Networks 23 (2010), 20-34
%
%     Watanabe(2010). Asymptotic Equivalance of Bayes Cross Validation and
%     Widely applicable Information Criterion in Singular Learning Theory.
%     Journal of Machine Learning Research 11 (2010), 3571-3594.
%     
%

% Copyright (c) 2011-2013 Ville Tolvanen

  ip=inputParser;
  ip.FunctionName = 'GP_WAIC';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('method', 'V', @(x) ismember(x,{'V' 'G'}))
  ip.addParamValue('form', 'mean', @(x) ismember(x,{'mean','all'}))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  method=ip.Results.method;
  form=ip.Results.form;
  % pass these forward
  options=struct();
  z = ip.Results.z;
  if ~isempty(ip.Results.z)
    options.zt=ip.Results.z;
    options.z=ip.Results.z;
  end
  
  [tn, nin] = size(x);
  
  % ====================================================
  if isstruct(gp)     % Single GP or MCMC solution
    switch gp.type
      case {'FULL' 'VAR' 'DTC' 'SOR'}
        tstind = [];
      case {'FIC' 'CS+FIC'}
        tstind = 1:tn;
      case 'PIC'
        tstind = gp.tr_index;
    end

    if isfield(gp, 'etr')
      % MCMC solution
      [Ef, Varf, BUt] = gpmc_preds(gp,x,y, x, 'yt', y, 'tstind', tstind, options);
      BUt=log(mean(exp(BUt),2));
      GUt = zeros(tn,1);
      Elog = zeros(tn,1);
      Elog2 = zeros(tn,1);
      
      nsamples = length(gp.edata);
      if strcmp(gp.type, 'PIC')
        tr_index = gp.tr_index;
        gp = rmfield(gp, 'tr_index');
      else
        tr_index = [];
      end
      
      %Ef = zeros(tn, nsamples);
      %Varf = zeros(tn, nsamples);
      sigma2 = zeros(tn, nsamples);
      for j = 1:nsamples
        Gp = take_nth(gp,j);
        if  strcmp(gp.type, 'FIC') | strcmp(gp.type, 'PIC')  || strcmp(gp.type, 'CS+FIC') || strcmp(gp.type, 'VAR') || strcmp(gp.type, 'DTC') || strcmp(gp.type, 'SOR')
          Gp.X_u = reshape(Gp.X_u,length(Gp.X_u)/nin,nin);
        end
        Gp.tr_index = tr_index;

        gp_array{j} = Gp;
        %[Ef(:,j), Varf(:,j)] = gp_pred(Gp, x, y, x, 'yt', y, 'tstind', tstind, options);
        if isfield(gp.lik.fh,'trcov')
          sigma2(:,j) = repmat(Gp.lik.sigma2,1,tn);
        end
      end
      
      if isequal(method,'V')
        % Evaluate WAIC using the Variance method
        
        if isfield(gp.lik.fh,'trcov')
          % Gaussian likelihood
          for i=1:tn
%             fmin = mean(Ef(i,:) - 9*sqrt(Varf(i,:)));
%             fmax = mean(Ef(i,:) + 9*sqrt(Varf(i,:)));
%             Elog(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
%                                        .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))').^2), fmin, fmax);
%             Elog2(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
%                                         .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))')), fmin, fmax);
%                                       
            m = Ef(i,:);
            s2 = Varf(i,:);
            m0 = 1; m1 = m; m2 = m.^2 + s2; m3 = m.*(m.^2+3*s2);
            m4 = m.^4+6.*m.^2.*s2+3*s2.^2;
            Elog2(i) = mean((-0.5.*log(2.*pi.*sigma2(i,:)) - y(i).^2./(2.*sigma2(i,:))).*m0 - 1./(2.*sigma2(i,:)) .* m2 + y(i)./sigma2(i,:) .* m1);
            Elog(i) = mean((1/4 .* m4 - y(i) .* m3 + (3.*y(i).^2./2+0.5.*log(2.*pi.*sigma2(i,:)).*sigma2(i,:)) .* m2 ...
              - (y(i).^3 + y(i).*log(2.*pi.*sigma2(i,:)).*sigma2(i,:)) .* m1 + (y(i).^4./4 + 0.5.*y(i).^2.*log(2.*pi.*sigma2(i,:)).*sigma2(i,:) ...
              + 0.25.*log(2.*pi.*sigma2(i,:)).^2.*sigma2(i,:).^2) .* m0) ./ sigma2(i,:).^2);
          end
          Elog2 = Elog2.^2;
          Vn = (Elog-Elog2);
          if strcmp(form, 'mean')
            Vn = mean(Vn);
            BUt = mean(BUt);
          end
          waic = BUt - Vn;
        else
          % non-Gaussian likelihood
          for i=1:tn
            if ~isempty(z)
              z1 = z(i);
            else
              z1 = [];
            end
            if ~isequal(gp.lik.type, 'Coxph')
              fmin = mean(Ef(i,:) - 9*sqrt(Varf(i,:)));
              fmax = mean(Ef(i,:) + 9*sqrt(Varf(i,:)));
              Elog(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
                .*llvec(gp_array, y(i), f, z1).^2), fmin, fmax);
              Elog2(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
                .*llvec(gp_array, y(i), f, z1)), fmin, fmax);
            else
              ntime = size(gp.lik.xtime,1);
              for i2=1:nsamples
                % Use MC to integrate over latents
                ns = 10000;
                Sigma_tmp = diag(Varf([1:ntime ntime+i],i2));
                f = mvnrnd(Ef([1:ntime ntime+i],i2), Sigma_tmp, ns);
                tmp2(i2) =  1/ns * sum(llvec(gp_array{i2}, y(i,:), f', z1));
                tmp(i2) = 1/ns * sum((llvec(gp_array{i2}, y(i,:), f', z1)).^2);
              end
              Elog2(i)=mean(tmp2);
              Elog(i)=mean(tmp);
            end
          end
          Elog2 = Elog2.^2;
          Vn = (Elog-Elog2);
          if strcmp(form, 'mean')
            Vn = mean(Vn);
            BUt = mean(BUt);
          end
          waic = BUt - Vn;
        end
        
      else
        % Evaluate WAIC using the expected value form via Gibbs training
        % loss
        
        if isfield(gp.lik.fh,'trcov')
          % Gaussian likelihood
          for i=1:tn
            fmin = mean(Ef(i,:) - 9*sqrt(Varf(i,:)));
            fmax = mean(Ef(i,:) + 9*sqrt(Varf(i,:)));
            GUt(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
                                      .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))')), fmin, fmax);
          end
          if strcmp(form, 'mean')
            GUt = mean(GUt);
            BUt = mean(BUt);
          end
          waic = BUt-2*(BUt-GUt);
        else
          % non-Gaussian likelihood
          for i=1:tn
            if ~isempty(z)
              z1 = z(i);
            else
              z1 = [];
            end
            fmin = mean(Ef(i,:) - 9*sqrt(Varf(i,:)));
            fmax = mean(Ef(i,:) + 9*sqrt(Varf(i,:)));
            GUt(i) = quadgk(@(f) mean(multi_npdf(f,Ef(i,:),(Varf(i,:))) ...
                                      .*llvec(gp_array, y(i), f, z1)), fmin, fmax);
          end
          if strcmp(form, 'mean')
            GUt = mean(GUt);
            BUt = mean(BUt);
          end
          waic = BUt-2*(BUt-GUt);
        end
      end
      
      
    else
      % A single GP solution
      [Ef, Varf, BUt] = gp_pred(gp, x, y, x, 'yt', y, 'tstind', tstind, options);

      GUt = zeros(tn,1);
      Elog = zeros(tn,1);
      Elog2 = zeros(tn,1);

      if isequal(method,'V')
        % Estimate WAIC with variance form
        
        if isfield(gp.lik.fh,'trcov')
          % Gaussian likelihood
          sigma2 = gp.lik.sigma2;
          
          for i=1:tn
            
            % Analytical moments for Gaussian distribution
            
            m0 = 1; m1 = Ef(i); m2 = Ef(i)^2 + Varf(i); m3 = Ef(i)*(Ef(i)^2+3*Varf(i));
            m4 = Ef(i)^4+6*Ef(i)^2*Varf(i)+3*Varf(i)^2;
          
            Elog2(i) = (-0.5*log(2*pi*sigma2) - y(i).^2./(2.*sigma2))*m0 - 1./(2.*sigma2) * m2 + y(i)./sigma2 * m1;
            Elog(i) = (1/4 * m4 - y(i) * m3 + (3*y(i).^2./2+0.5*log(2*pi*sigma2).*sigma2) * m2 ...
                       - (y(i).^3 + y(i).*log(2*pi*sigma2).*sigma2) * m1 + (y(i).^4/4 + 0.5*y(i).^2*log(2*pi*sigma2).*sigma2 ...
                                                              + 0.25*log(2*pi*sigma2).^2.*sigma2.^2) * m0) ./ sigma2.^2;
            
          end
          Elog2 = Elog2.^2;
          Vn = Elog-Elog2;
          if strcmp(form,'mean')
            BUt = mean(BUt);
            Vn = mean(Vn);
          end
          waic = BUt - Vn;

        else
          % Non-Gaussian likelihood
          for i=1:tn
            if ~isempty(z)
              z1 = z(i);
            else
              z1 = [];
            end
            if ~isequal(gp.lik.type, 'Coxph')
              fmin = Ef(i)-9*sqrt(Varf(i));
              fmax = Ef(i)+9*sqrt(Varf(i));
              Elog(i) = quadgk(@(f) norm_pdf(f, Ef(i), sqrt(Varf(i))).*llvec(gp, y(i), f, z1).^2 ,...
                fmin, fmax);
              Elog2(i) = quadgk(@(f) norm_pdf(f, Ef(i), sqrt(Varf(i))).*llvec(gp, y(i), f, z1) ,...
                fmin, fmax);
            else
              % Use MC to integrate over latents
              ntime = size(gp.lik.xtime,1);
              ns = 10000;
              Sigma_tmp = Varf([1:ntime ntime+i], [1:ntime ntime+i]);
              Sigma_tmp = (Sigma_tmp + Sigma_tmp') ./ 2;
              f = mvnrnd(Ef([1:ntime ntime+i]), Sigma_tmp, ns);
              Elog2(i) = 1/ns * sum(llvec(gp, y(i,:), f', z1));
              Elog(i) = 1/ns * sum((llvec(gp, y(i,:), f', z1)).^2);
            end
          end
          Elog2 = Elog2.^2;
          Vn = Elog-Elog2;
          if strcmp(form, 'mean')
            Vn = mean(Vn);
            BUt = mean(BUt);
          end
          waic = BUt - Vn;
          
        end
        
      else
        % WAIC using the expected value form via Gibbs training loss GUt
        
        if isfield(gp.lik.fh,'trcov')
          % Gaussian likelihood
          sigma2 = gp.lik.sigma2;
          for i=1:tn
            if Varf(i)<eps
              GUt(i)=(-0.5*log(2*pi*sigma2)- (y(i) - Ef(i)).^2/(2.*sigma2));
            else
              
              % GUt(i) = quadgk(@(f) norm_pdf(f,Ef(i),sqrt(Varf(i))).*(-0.5*log(2*pi*sigma2)- (y(i) - f).^2/(2.*sigma2)), fmin, fmax);

              m0 = 1; m1 = Ef(i); m2 = Ef(i)^2 + Varf(i);
              
              GUt(i) = (-0.5*log(2*pi*sigma2) - y(i).^2./(2.*sigma2))*m0 - 1./(2.*sigma2) * m2 + y(i)./sigma2 * m1;
            end
          end
          if strcmp(form,'mean')
            GUt = mean(GUt);
            BUt = mean(BUt);
          end
          waic = BUt-2*(BUt-GUt);
        else
          % Non-Gaussian likelihood
          for i=1:tn
            if ~isempty(z)
              z1 = z(i);
            else
              z1 = [];
            end
            if ~isequal(gp.lik.type, 'Coxph')
              fmin = Ef(i)-9*sqrt(Varf(i));
              fmax = Ef(i)+9*sqrt(Varf(i));
              GUt(i) = quadgk(@(f) norm_pdf(f, Ef(i), sqrt(Varf(i))).*llvec(gp, y(i), f, z1) ,...
                fmin, fmax);
            else
              % If likelihood coxph use mc to integrate over latents
              ntime = size(gp.lik.xtime,1);
              ns = 10000;
              Sigma_tmp = Varf([1:ntime ntime+i], [1:ntime ntime+i]);
              Sigma_tmp = (Sigma_tmp + Sigma_tmp') ./ 2;
              f = mvnrnd(Ef([1:ntime ntime+i]), Sigma_tmp, ns);
              GUt(i) = 1/ns * sum(llvec(gp, y(i), f', z1));
            end
          end
          if strcmp(form,'mean')
            GUt = mean(GUt);
            BUt = mean(BUt);
          end
          waic = BUt-2*(BUt-GUt);
        end
        
      end
      
      
    end
    
  elseif iscell(gp)
    
    % gp_ia solution
    
    switch gp{1}.type
      case {'FULL' 'VAR' 'DTC' 'SOR'}
        tstind = [];
      case {'FIC' 'CS+FIC'}
        tstind = 1:tn;
      case 'PIC'
        tstind = gp{1}.tr_index;
    end
    
    [tmp, tmp, BUt] = gp_pred(gp,x,y, x, 'yt', y, 'tstind', tstind, options);
    GUt = zeros(tn,1);
    Elog = zeros(tn,1);
    Elog2 = zeros(tn,1);
    
    nsamples = length(gp);
    for j = 1:nsamples
      Gp = gp{j};
      weight(j) = Gp.ia_weight;
      w(j,:) = gp_pak(Gp);
      [Ef(:,j), Varf(:,j)] = gp_pred(Gp, x, y, x, 'yt', y, 'tstind', tstind, options);
      if isfield(Gp.lik.fh,'trcov')
        sigma2(:,j) = repmat(Gp.lik.sigma2,1,tn);
      end
    end
    if isequal(method,'V')
      % Evaluate WAIC using the variance form
      
      if isfield(gp{1}.lik.fh,'trcov')
        % Gaussian likelihood
        for i=1:tn
          fmin = sum(weight.*Ef(i,:) - 9*weight.*sqrt(Varf(i,:)));
          fmax = sum(weight.*Ef(i,:) + 9*weight.*sqrt(Varf(i,:)));
          Elog(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                    .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))').^2), fmin, fmax);
          Elog2(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                     .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))')), fmin, fmax);
        end
        Elog2 = Elog2.^2;
        Vn = (Elog-Elog2);
        if strcmp(form, 'mean')
          Vn = mean(Vn);
          BUt = mean(BUt);
        end
        waic = BUt - Vn;
      else
        % non-Gaussian likelihood
        for i=1:tn
          if ~isempty(z)
            z1 = z(i);
          else
            z1 = [];
          end
          fmin = sum(weight.*Ef(i,:) - 9*weight.*sqrt(Varf(i,:)));
          fmax = sum(weight.*Ef(i,:) + 9*weight.*sqrt(Varf(i,:)));
          Elog(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                    .*llvec(gp, y(i), f, z1).^2), fmin, fmax);
          Elog2(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                     .*llvec(gp, y(i), f, z1)), fmin, fmax);
        end
        Elog2 = Elog2.^2;
        Vn = (Elog-Elog2);
        if strcmp(form, 'mean')
          Vn = mean(Vn);
          BUt = mean(BUt);
        end
        waic = BUt - Vn;
        
      end
      
    else
      % Evaluate WAIC using the expected value form via Gibbs training loss
      
      if isfield(gp{1}.lik.fh,'trcov')
        % Gaussian likelihood
        for i=1:tn
          fmin = sum(weight.*Ef(i,:) - 9*weight.*sqrt(Varf(i,:)));
          fmax = sum(weight.*Ef(i,:) + 9*weight.*sqrt(Varf(i,:)));
          GUt(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                   .*bsxfun(@minus,-bsxfun(@rdivide,(repmat((y(i)-f),nsamples,1)).^2,(2.*sigma2(i,:))'), 0.5*log(2*pi*sigma2(i,:))')), fmin, fmax);
        end
        if strcmp(form, 'mean')
          GUt = mean(GUt);
          BUt = mean(BUt);
        end
        waic = BUt-2*(BUt-GUt);

      else
        % non-gaussian likelihood
        for i=1:tn
          if ~isempty(z)
            z1 = z(i);
          else
            z1 = [];
          end
          fmin = sum(weight.*Ef(i,:) - 9*weight.*sqrt(Varf(i,:)));
          fmax = sum(weight.*Ef(i,:) + 9*weight.*sqrt(Varf(i,:)));
          GUt(i) = quadgk(@(f) sum(bsxfun(@times, multi_npdf(f,Ef(i,:),(Varf(i,:))),weight') ...
                                   .*llvec(gp, y(i), f, z1)), fmin, fmax);
        end
        if strcmp(form, 'mean')
          GUt = mean(GUt);
          BUt = mean(BUt);
        end
        waic = BUt-2*(BUt-GUt);

      end
    end
    
  end

end

function lls=llvec(gp, y, fs, z)
% Compute a vector of lls for vector argument fs used by quadgk. In
% case of IA or MC, return a matrix with rows corresponding to one
% GP and columns corresponding to all of the GP's.
  
  if isstruct(gp)
      % single gp
      lls=zeros(1,size(fs,2));
      for i1=1:size(fs,2)
        lls(i1)=gp.lik.fh.ll(gp.lik,y,fs(:,i1),z);
      end
      %     else
      %       % mc
      %       lls=zeros(length(gp), length(fs));
      %       for i=1:numel(fs)
      %         for j=1:numel(gp.edata)
      %           Gp = take_nth(gp, j);
      %           lls(j,i) = Gp.lik.fh.ll(Gp.lik, y, fs(i), z);
      %         end
      %       end
  else
    % ia & mc
    lls=zeros(length(gp), length(fs));
    for i=1:numel(fs)
      for j=1:numel(gp)
        lls(j,i) = gp{j}.lik.fh.ll(gp{j}.lik, y, fs(i), z);
      end  
    end
  end
end

function mpdf = multi_npdf(f, mean, sigma2)
% for every element in f, compute means calculated with 
% norm_pdf(f(i), mean, sqrt(sigma2)). If mean and sigma2
% are vectors, returns length(mean) x length(f) matrix. 
  
  mpdf = zeros(length(mean), length(f));
  for i=1:length(f)
    mpdf(:,i) = norm_pdf(f(i), mean, sqrt(sigma2));
  end
end

function [m_0, m_1, m_2, m_3, m_4] = moments(fun, a, b, rtol, atol, minsubs)
% QUAD_MOMENTS Calculate the 0th, 1st and 2nd moment of a given
%              (unnormalized) probability distribution
%
%   [m_0, m_1, m_2] = quad_moments(fun, a, b, varargin) 
%   Inputs:
%      fun  = Function handle to the unnormalized probability distribution
%      a,b  = integration limits [a,b]
%      rtol = relative tolerance for the integration (optional, default 1e-6)
%      atol = absolute tolerance for the integration (optional, default 1e-10)
%               
%   Returns the first three moments:
%      m0  = int_a^b fun(x) dx
%      m1  = int_a^b x*fun(x) dx / m0
%      m2  = int_a^b x^2*fun(x) dx / m0
%
%   The function uses an adaptive Gauss-Kronrod quadrature. The same set of 
%   integration points and intervals are used for each moment. This speeds up 
%   the evaluations by factor 3, since the function evaluations are done only 
%   once.
% 
%   The quadrature method is described by:
%   L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
%   Journal of Computational and Applied Mathematics, 211, 2008, 
%   pp. 131-140.

%   Copyright (c) 2010 Jarno Vanhatalo, Jouni Hartikainen
  
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  maxsubs = 650;
  
  if nargin < 4
    rtol = 1.e-6;
  end
  if nargin < 5
    atol = 1.e-10;
  end
  if nargin < 6
    minsubs = 10;
  end
  
  rtol = max(rtol,100*eps);
  atol = max(atol,0);
  minsubs = max(minsubs,2); % At least two subintervals are needed
  
  % points and weights
  points15 = [0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
              0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
              0.9914553711208126];
  points = [-points15(end:-1:1); 0; points15];
  
  w15 = [0.2044329400752989, 0.1903505780647854, 0.1690047266392679, ...
         0.1406532597155259, 0.1047900103222502, 0.06309209262997855, ...
         0.02293532201052922];
  w = [w15(end:-1:1), 0.2094821410847278, w15];
  
  w7 = [0,0.3818300505051189,0,0.2797053914892767,0,0.1294849661688697,0];
  ew = w - [w7(end:-1:1), 0.4179591836734694, w7];
  
  samples = numel(w);
  
  % split the interval.
  if b-a <= 0
    c = a; a = b; b=c;
    warning('The start of the integration interval was less than the end of it.')
  end
  apu = a + (1:(minsubs-1))./minsubs*(b-a);
  apu = [a,apu,b];
  subs = [apu(1:end-1);apu(2:end)];
  
  % Initialize partial sums.
  Ifx_ok = 0;
  Ifx1_ok = 0;
  Ifx2_ok = 0;
  Ifx3_ok = 0;
  Ifx4_ok = 0;
  % The main loop
  while true
    % subintervals and their midpoints
    midpoints = sum(subs)/2;   
    halfh = diff(subs)/2;  
    x = bsxfun(@plus,points*halfh,midpoints);
    x = reshape(x,1,[]);
    
    fx = fun(x);
    fx1 = fx.*x;
    fx2 = fx.*x.^2;
    fx3 = fx.*x.^3;
    fx4 = fx.*x.^4;

    fx = reshape(fx,samples,[]);
    fx1 = reshape(fx1,samples,[]);
    fx2 = reshape(fx2,samples,[]);
    fx3 = reshape(fx3,samples,[]);
    fx4 = reshape(fx4,samples,[]);
    
    % Subintegrals.
    Ifxsubs = (w*fx) .* halfh;
    errsubs = (ew*fx) .* halfh;
    Ifxsubs1 = (w*fx1) .* halfh;
    Ifxsubs2 = (w*fx2) .* halfh;
    Ifxsubs3 = (w*fx3) .* halfh;
    Ifxsubs4 = (w*fx4) .* halfh;

    % Ifx and tol.
    Ifx = sum(Ifxsubs) + Ifx_ok;
    Ifx1 = sum(Ifxsubs1) + Ifx1_ok;
    Ifx2 = sum(Ifxsubs2) + Ifx2_ok;
    Ifx3 = sum(Ifxsubs3) + Ifx3_ok;
    Ifx4 = sum(Ifxsubs4) + Ifx4_ok;
    tol = max(atol,rtol*abs(Ifx));
    
    % determine the indices ndx of Ifxsubs for which the
    % errors are acceptable and remove those from subs
    ndx = find(abs(errsubs) <= (2/(b-a)*halfh*tol));
    subs(:,ndx) = [];
    if isempty(subs)
      break
    end
    
    % Update the integral.
    Ifx_ok = Ifx_ok + sum(Ifxsubs(ndx));
    Ifx1_ok = Ifx1_ok + sum(Ifxsubs1(ndx));
    Ifx2_ok = Ifx2_ok + sum(Ifxsubs2(ndx));
    Ifx3_ok = Ifx3_ok + sum(Ifxsubs3(ndx));
    Ifx4_ok = Ifx4_ok + sum(Ifxsubs4(ndx));

    
    % Quit if too many subintervals.
    nsubs = 2*size(subs,2);
    if nsubs > maxsubs
      warning('quad_moments: Reached the limit on the maximum number of intervals in use.');
      break
    end
    midpoints(ndx) = []; 
    subs = reshape([subs(1,:); midpoints; midpoints; subs(2,:)],2,[]); % Divide the remaining subintervals in half
  end
  
  % Scale moments
  m_0 = Ifx;
  m_1 = Ifx1./Ifx;
  m_2 = Ifx2./Ifx;
  m_3 = Ifx3./Ifx;
  m_4 = Ifx4./Ifx;
end
