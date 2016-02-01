function [dic, p_eff, Davg] = gp_dic(gp, x, y, varargin);
%GP_DIC The DIC and the effective number of parameters in a GP model
%
%  Description
%    [DIC, P_EFF] = GP_DIC(GP, X, Y) evaluates DIC and the effective
%    number of parameters as defined by Spiegelhalter et al
%    (2002). The statistics are evaluated with focus on parameters
%    or latent variables depending on the input GP (See
%    Spiegelhalter et al (2002) for discussion on the parameters
%    in focus in Bayesian model). X contains training inputs and Y
%    training outputs. Output form can be set by providing 
%    'output', PARAM ('DIC', 'mlpd') and 'form', PARAM ('mean', 'all').
%
%    DIC and p_eff are evaluated as follows:
%     1) GP is a record structure from gp_mc or an array of GPs
%        from gp_ia. In this case the focus is in the parameters
%        (the parameters of the covariance function and the
%        likelihood). The DIC and the effective number of
%        parameters are evaluated as described in equation (6.10)
%        of Bayesian Data Analysis, second edition (Gelman et al):
%          p_eff = E[D(y, th)|y] - D(y, E[th|y])
%          DIC   = E[D(y, th)|y] + p_eff
%        where all the expectations are taken over p(th|y) and 
%          D(y, th) = -2log(p(y|th)). 
%        Now in this formulation we first marginalize over the
%        latent variables to obtain 
%          p(y|th) = \int p(y|f)(p(f|th) df. 
%        If the likelihood is non-Gaussian the marginalization can
%        not be performed exactly. In this case, if GP is an MCMC
%        record, we use Laplace approximation to approximate
%        p(y|th). If GP is IA array, we use the either EP or
%        Laplace approximation, depending which has been used in
%        gp_ia.
%
%     2) GP is Gaussian process structure
%
%        In this case the focus is in the latent variables and the
%        parameters are considered fixed. The mean of the
%        deviance is now evaluated as
%               E[D(y, f)|y] = -2 \int log[p(y|f)] p(f|th) df
%
%     3) GP is a record structure from gp_mc or an array of GPs from 
%        gp_ia, but the focus is defined to be both latent-variables 
%        and parameters, 
%               [DIC, P_EFF] = GP_DIC(GP, X, Y, 'focus', 'all')
%        In this case the focus will be the latent variables and
%        parameters. Thus now we will use the posterior p(f, th|y)
%        instead of the conditional posterior p(f|th,y) or
%        posterior marginal p(th|y). The DIC and the effective
%        number of parameters are evaluated as described in
%        equation (6.10) of Bayesian Data Analysis, second edition
%        (Gelman et al):
%               p_eff = E[D(y, f, th)|y] - D(y, E[f, th|y])
%               DIC   = E[D(y, f, th)|y] + p_eff
%        where all the expectations are taken over p(f,th|y).
%       
%  See also
%    GP_PEFF, DEMO_MODELASSESMENT1
%
%  References: 
%
%    Spiegelhalter, Best, Carlin and van der Linde (2002). Bayesian
%    measures of model complexity and fit. J. R. Statist. Soc. B,
%    64, 583-639.
%         
%    Gelman, Carlin, Stern and Rubin (2004) Bayesian Data Analysis,
%    second edition. Chapman & Hall / CRC.
%   

% Copyright (c) 2009-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.xt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_DIC';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('focus', 'param', @(x) ismember(x,{'param','latent','all'}))
  ip.addParamValue('output', 'DIC', @(x) ismember(x,{'DIC', 'mlpd'}))
  ip.addParamValue('form', 'mean', @(x) ismember(x,{'mean','all'}))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  focus=ip.Results.focus;
  output=ip.Results.output;
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
      if nargin < 4 || isempty(focus)
        focus = 'param';
      end
    else
      % A single GP
      focus = 'latent';
    end     
    
    fh_pred = @gp_pred;
    
    switch focus
      
      case 'param'
        % An MCMC solution and focus in parameters
        
        % evaluate the mean of the parameters
        if strcmp(gp.type, 'PIC')
          tr_index = gp.tr_index;
          gp = rmfield(gp, 'tr_index');
        else
          tr_index = [];
        end
        for i = 1:length(gp.edata)
          Gp = take_nth(gp,i);
          w(i,:) = gp_pak(Gp);
        end                
        Gp.tr_index = tr_index;
        Gp = gp_unpak(Gp, mean(w,1));
        if strcmp(gp.type, 'FIC') | strcmp(gp.type, 'PIC')  || strcmp(gp.type, 'CS+FIC') || strcmp(gp.type, 'VAR') || strcmp(gp.type, 'DTC') || strcmp(gp.type, 'SOR')
          Gp.X_u = reshape(Gp.X_u,length(Gp.X_u)/nin,nin);
        end

        if isfield(gp.lik.fh,'trcov')
          % a Gaussian likelihood
          if strcmp(form, 'mean')
            Davg = 2*mean(gp.edata);
            [e, edata] = gp_e(mean(w,1), Gp, x, y);
            Dth = 2*edata;
          else
            Davg = 2.*gp.edata;
            Dth = zeros(size(Davg));
            for i1=1:size(w,1)
              [e, edata] = gp_e(w(i,:), Gp, x, y);
              Dth(i1) = 2*edata;
            end
          end
        else
          % non-Gaussian likelihood
          
          % For non-Gaussian likelihood we cannot evaluate the marginal
          % likelihood p(y|th) exactly. For this reason we
          % use Laplace approximation to approximate p(y|th)
          gp2 = Gp;
          gp2 = gp_set(gp2, 'latent_method', 'Laplace');
          if strcmp(form, 'mean')
            [e, edata] = gpla_e(mean(w,1), gp2, x, y, 'z', z);
            Dth = 2.*edata;
            
            for i1 = 1:length(gp.edata)
              [e, edata] = gpla_e(w(i1,:), gp2, x, y, 'z', z);
              Davg(i1) = 2.*edata;
            end
            Davg = mean(Davg);
          else
            Dth = zeros(size(gp.edata,1));
            Davg = zeros(size(gp.edata,1));
            for i1 = 1:length(gp.edata)
              [tmp, edata] = gpla_e(w(i1,:), gp2, x, y, 'z', z);
              Dth(i1) = 2.*edata;
              [tmp, edata] = gpla_e(w(i1,:), gp2, x, y, 'z', z);
              Davg(i1) = 2.*edata;
            end
          end
        end
        
      case 'latent'     
        % A single GP solution -> focus on latent variables

        if ~isfield(gp.lik, 'type_nd')
          [Ef, Varf, lpy, Ey, VarY] = fh_pred(gp, x, y, x, 'yt', y, 'tstind', tstind, options);
          sampf = gp_rnd(gp, x, y, x, 'tstind', tstind, 'nsamp', 5000, options);
        else
          [Ef, Covf, lpy, Ey, Covy] = gp_jpred(gp, x, y, x, 'yt', y, 'tstind', tstind, options);
          sigma_tmp = chol((Covf+Covf')./2, 'lower');
          sampf = repmat(Ef,1,5000) + sigma_tmp*randn(size(Ef,1),5000);          
          Varf = diag(Covf);
          VarY = diag(Covy);
        end
        if isfield(gp.lik.fh,'trcov')
          % a Gaussian likelihood
          sigma2 = VarY - Varf;
          if strcmp(form, 'mean')
            Dth = sum(log(2*pi*sigma2)) + sum( (y - Ef).^2./sigma2 );
            Davg = sum(log(2*pi*sigma2)) + mean(sum( (repmat(y,1,5000) - sampf).^2./repmat(sigma2,1,5000), 1));
          else
            Dth = log(2*pi*sigma2) + (y - Ef).^2./sigma2;
            Davg = log(2*pi*sigma2) + mean((repmat(y,1,5000) - sampf).^2./repmat(sigma2,1,5000),2);  
          end
        else 
          % non-Gaussian likelihood
          if strcmp(form, 'mean')
            Dth = -2.*gp.lik.fh.ll(gp.lik, y, Ef, z);
            for i1 = 1:size(sampf, 2)
              Davg(i1) = gp.lik.fh.ll(gp.lik, y, sampf(:,i1), z);
            end
            Davg = -2.*mean(Davg);
          else
            if ~isequal(gp.lik.type, 'Coxph')
              if isempty(z)
                z = ones(size(Ef));
              end
              Dth = -2.*arrayfun(@(a,b,c) gp.lik.fh.ll(gp.lik, a, b, c), y, Ef, z);
              for i1 = 1:size(sampf, 2)
                Davg(:,i1) = arrayfun(@(a,b,c) gp.lik.fh.ll(gp.lik, a, b, c), y, sampf(:,i1), z);
              end
              Davg = -2.*mean(Davg,2);
            else
              % If likelihood coxph use mc to integrate over latents
              ntime = size(gp.lik.xtime,1);
              for i=1:tn
                Dth(i,1) = -2*gp.lik.fh.ll(gp.lik, y(i), Ef([1:ntime ntime+i]), z(i));
                fs = repmat(Ef([1:ntime ntime+i]),1,5000) + chol((Covf([1:ntime ntime+i],[1:ntime ntime+i])+Covf([1:ntime ntime+i],[1:ntime ntime+i])')./2, 'lower')*randn(ntime+1,5000);
                for i1=1:size(fs,1)
                  Davg(i,i1) = gp.lik.fh.ll(gp.lik, y(i), fs(:,i1), z(i));
                end
              end
              Davg = -2.*mean(Davg,2);
            end
          end
        end
        
      case 'all'        
        % An MCMC solution and focus on all parameters
        
        nsamples = length(gp.edata);
        if strcmp(gp.type, 'PIC')
          tr_index = gp.tr_index;
          gp = rmfield(gp, 'tr_index');
        else
          tr_index = [];
        end
        for i = 1:nsamples
          Gp = take_nth(gp,i);
          w(i,:) = gp_pak(Gp);
          if  strcmp(gp.type, 'FIC') | strcmp(gp.type, 'PIC')  || strcmp(gp.type, 'CS+FIC') || strcmp(gp.type, 'VAR') || strcmp(gp.type, 'DTC') || strcmp(gp.type, 'SOR')
            Gp.X_u = reshape(Gp.X_u,length(Gp.X_u)/nin,nin);
          end
          Gp.tr_index = tr_index;
          if isfield(gp.lik.fh,'trcov')
            % a Gaussian likelihood
            sampf(:,i) = gp_rnd(Gp, x, y, x, 'tstind', tstind, options);
            [Ef(:,i), Varf, lpy, Ey, VarY] = gp_pred(Gp, x, y, x, 'yt', y, 'tstind', tstind);
            sigma2(:,i) = VarY - Varf;
          end
        end
        
        if isfield(gp.lik.fh,'trcov')
          % a Gaussian likelihood
          if strcmp(form, 'mean')
            Ef = mean(Ef, 2);
            msigma2 = mean(sigma2,2);
            Dth = sum(log(2*pi*msigma2)) + sum( (y - Ef).^2./msigma2 );
            Davg = mean(sum(log(2*pi*sigma2),1)) + mean(sum( (repmat(y,1,nsamples) - sampf).^2./sigma2, 1));
          else
            Ef = mean(Ef, 2);
            msigma2 = mean(sigma2,2);
            Dth = (log(2*pi*msigma2)) + ( (y - Ef).^2./msigma2 );
            Davg = mean(log(2*pi*sigma2),2) + mean( (repmat(y,1,nsamples) - sampf).^2./sigma2, 2);
          end
        else 
          % non-Gaussian likelihood
          Gp = gp_unpak(Gp, mean(w,1));
          if strcmp(form, 'mean')
            Dth = -2.*Gp.lik.fh.ll(Gp.lik, y, mean(gp.latentValues,1)', z);
            for i1 = 1:nsamples
              Gp = take_nth(gp,i1);
              Davg(i1) = Gp.lik.fh.ll(Gp.lik, y, Gp.latentValues', z);
            end
            Davg = -2.*mean(Davg);
          else
            if ~isempty(z)
              z1 = z;
            else
              z1 = ones(size(y));
            end
            Dth = -2.*arrayfun(@(a,b,c) Gp.lik.fh.ll(Gp.lik, a, b, c), y, mean(gp.latentValues,1)',z1);
            for i1 = 1:nsamples
              Gp = take_nth(gp,i1);
              Davg(:,i1) = arrayfun(@(a,b,c) Gp.lik.fh.ll(Gp.lik, a, b, c), y, Gp.latentValues', z1);
            end
            Davg = -2.*mean(Davg,2);
          end
        end
    end       

    % ====================================================
  elseif iscell(gp)
    % gp_ia solution
    if nargin < 4 || isempty(focus)
      focus = 'param';
    end
    if strcmp(focus, 'latent')
      error(['gp_dic: The focus can be ''latent'' only if single GP structure is given. '...
             'With IA cell array possible options are ''param'' and ''all''.            ']);
    end
    
    switch gp{1}.type
      case {'FULL' 'VAR' 'DTC' 'SOR'}
        tstind = [];
      case {'FIC' 'CS+FIC'}
        tstind = 1:tn;
      case 'PIC'
        tstind = gp{1}.tr_index;
    end

    % Define the error and prediction functions
    if isstruct(gp{1}.lik) && isfield(gp{1}, 'latent_method')
      fh_pred = gp{1}.fh.pred;
      fh_e = gp{1}.fh.e;
    else
      fh_e = @gp_e;
      fh_pred = @gp_pred;
    end
    
    switch focus
      
      case 'param'
        % An IA solution and focus in parameters
        
        for i = 1:length(gp)
          Gp = gp{i};
          weight(i) = Gp.ia_weight; 
          w(i,:) = gp_pak(Gp);
          [e, edata] = fh_e(w(i,:), Gp, x, y, options);
          energy(i) = edata;
        end
        Davg = 2*sum(energy.*weight);
        wh = sum(w.*repmat(weight',1,size(w,2)),1);
        Gp = gp_unpak(Gp, wh);
        [e, edata] = fh_e(wh, Gp, x, y, options);
        Dth = 2*edata;
        
      case 'all'
        % An IA solution and focus on all parameters
        
        nsamples = length(gp);
        for i = 1:nsamples
          Gp = gp{i};
          weight(i) = Gp.ia_weight;
          w(i,:) = gp_pak(Gp);
          [Ef(:,i), Varf(:,i), lpy, Ey, VarY] = fh_pred(Gp, x, y, x, 'yt', y, 'tstind', tstind, options);
          sigma2(:,i) = VarY - Varf(:,i);
        end
        mEf = sum(Ef.*repmat(weight, size(Ef,1), 1), 2);

        if isfield(gp{1}.lik.fh,'trcov')
          % a Gaussian likelihood
          msigma2 = sum(sigma2.*repmat(weight, size(Ef,1), 1), 2);
          if strcmp(form, 'mean')
            Dth = sum(log(2*pi*msigma2)) + sum( (y - mEf).^2./msigma2 );
            deviance = sum(log(2*pi*sigma2),1) + sum((Varf+Ef.^2-2.*repmat(y,1,nsamples).*Ef+repmat(y.^2,1,nsamples))./sigma2,1);
            Davg = sum(deviance.*weight);
          else
            Dth = log(2*pi*msigma2) + (y - mEf).^2./msigma2;
            deviance = log(2*pi*sigma2) + (Varf+Ef.^2-2.*repmat(y,1,nsamples).*Ef+repmat(y.^2,1,nsamples))./sigma2;
            Davg = sum(bsxfun(@times,deviance,weight),2);
          end
        else
          % non-Gaussian likelihood
          mw = sum(w.*repmat(weight', 1, size(w,2)), 1);
          Gp = gp_unpak(Gp, mw);
          if strcmp(form, 'mean')
            Dth = -2.*Gp.lik.fh.ll(Gp.lik, y, mEf, z);
            for i1 = 1:nsamples
              Gp = gp{i1};
              Davg(i1) = Gp.lik.fh.ll(Gp.lik, y, Ef(:,i), z);
            end
            Davg = -2.*sum(Davg.*weight);
          else
            if isempty(z)
              z = zeros(size(y));
            end
            Dth = -2.*arrayfun(@(a,b,c) Gp.lik.fh.ll(Gp.lik, a, b, c), y, mEf, z);
            for i1 = 1:nsamples
              Gp = gp{i1};
              Davg(:,i1) = arrayfun(@(a,b,c) Gp.lik.fh.ll(Gp.lik, a, b, c), y, Ef(:,i), z);
            end
            Davg = -2.*sum(bsxfun(@times, Davg, weight),2);
          end
        end
    end       

    % ====================================================
  else 
    error(['gp_dic: The first input must be a GP structure, a record structure';
           'from gp_mc or an array of GPs from gp_ia.                         ']) 
  end

  dic = 2*Davg - Dth;
  p_eff = Davg - Dth;
  if strcmp(output,'mlpd')
    % Scales output to correspond with waic, refpred, etc.
    if strcmp(form, 'mean')
      dic = dic/(-2*tn);
    else
      dic = dic/(-2);
    end
  end
  
end
