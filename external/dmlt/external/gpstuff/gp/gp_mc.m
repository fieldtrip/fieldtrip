function [record, gp, opt] = gp_mc(gp, x, y, varargin)
%GP_MC  Markov chain Monte Carlo sampling for Gaussian process models
%
%  Description
%    [RECORD, GP, OPT] = GP_MC(GP, X, Y, OPTIONS) Takes the Gaussian 
%    process structure GP, inputs X and outputs Y. Returns record 
%    structure RECORD with parameter samples, the Gaussian process GP
%    at current state of the sampler and an options structure OPT 
%    containing all the options in OPTIONS and information of the
%    current state of the sampler (e.g. the random number seed)
%
%    OPTIONS is optional parameter-value pair
%      z           - Optional observed quantity in triplet (x_i,y_i,z_i).
%                    Some likelihoods may use this. For example, in
%                    case of Poisson likelihood we have z_i=E_i,
%                    that is, expected value for ith case.
%      repeat      - Number of iterations between successive sample saves
%                    (that is every repeat'th sample is stored). Default 1.
%      nsamples    - Number of samples to be returned
%      display     - Defines if sampling information is printed, 1=yes, 0=no.
%                    Default 1. If >1, only every nth iteration is displayed.
%      hmc_opt     - Options structure for HMC sampler (see hmc2_opt). 
%                    When this is given the covariance function and
%                    likelihood parameters are sampled with hmc2
%                    (respecting infer_params option). If optional
%                    argument hmc_opt.nuts = 1, No-U-Turn HMC is used 
%                    instead. With NUTS, only mandatory parameter is
%                    number of adaptation steps hmc_opt.nadapt of step-size 
%                    parameter. For additional info, see hmc_nuts.
%      sls_opt     - Options structure for slice sampler (see sls_opt). 
%                    When this is given the covariance function and
%                    likelihood parameters are sampled with sls
%                    (respecting infer_params option).
%      latent_opt  - Options structure for latent variable sampler. When this 
%                    is given the latent variables are sampled with
%                    function stored in the gp.fh.mc field in the
%                    GP structure. See gp_set.
%      lik_hmc_opt - Options structure for HMC sampler (see hmc2_opt). 
%                    When this is given the parameters of the
%                    likelihood are sampled with hmc2. This can be
%                    used to have different hmc options for
%                    covariance and likelihood parameters.
%      lik_sls_opt - Options structure for slice sampler (see sls_opt). 
%                    When this is given the parameters of the
%                    likelihood are sampled with hmc2. This can be
%                    used to have different hmc options for
%                    covariance and likelihood parameters.
%      lik_gibbs_opt
%                  - Options structure for Gibbs sampler. Some likelihood
%                    function parameters need to be sampled with
%                    Gibbs sampling (such as lik_smt). The Gibbs
%                    sampler is implemented in the respective lik_*
%                    file.
%      persistence_reset 
%                  - Reset the momentum parameter in HMC sampler after 
%                    every repeat'th iteration, default 0.  
%      record      - An old record structure from where the sampling is 
%                    continued
%         
%    The GP_MC function makes nsamples*repeat iterations and stores
%    every repeat'th sample. At each iteration it samples first the
%    latent variables (if 'latent_opt' option is given), then the
%    covariance and likelihood parameters (if 'hmc_opt', 'sls_opt'
%    or 'gibbs_opt' option is given and respecting infer_params
%    option), and for last the the likelihood parameters (if
%    'lik_hmc_opt' or 'lik_sls_opt' option is given).
%
%  See also:
%    DEMO_CLASSIFIC1, DEMO_ROBUSTREGRESSION

% Copyright (c) 1998-2000,2010 Aki Vehtari
% Copyright (c) 2007-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

%#function gp_e gp_g

  ip=inputParser;
  ip.FunctionName = 'GP_MC';
  ip.addRequired('gp',@isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('nsamples', 1, @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('repeat', 1, @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('display', 1, @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('record',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('hmc_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('sls_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('ssls_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('latent_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('lik_hmc_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('lik_sls_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('lik_gibbs_opt', [], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('persistence_reset', 0, @(x) ~isempty(x) && isreal(x));
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;
  opt.nsamples=ip.Results.nsamples;
  opt.repeat=ip.Results.repeat;
  opt.display=ip.Results.display;
  record=ip.Results.record;
  opt.hmc_opt = ip.Results.hmc_opt;
  opt.ssls_opt = ip.Results.ssls_opt;
  opt.sls_opt = ip.Results.sls_opt;
  opt.latent_opt = ip.Results.latent_opt;
  opt.lik_hmc_opt = ip.Results.lik_hmc_opt;
  opt.lik_sls_opt = ip.Results.lik_sls_opt;
  opt.lik_gibbs_opt = ip.Results.lik_gibbs_opt;
  opt.persistence_reset = ip.Results.persistence_reset;
  
%   if isfield(gp.lik, 'nondiagW');
%     switch gp.lik.type
%       case {'LGP', 'LGPC'}
%         error('gp2_mc not implemented for this type of likelihood');
%       case {'Softmax', 'Multinom'}
%         [n,nout] = size(y);
%       otherwise
%         n = size(y,1);
%         nout=length(gp.comp_cf);
%     end
%   end
  
  
  % Default samplers and some checking
  if isfield(gp,'latent_method') && isequal(gp.latent_method,'MCMC')
    % If no options structures, use SSLS as a default sampler for hyperparameters
    % and ESLS for latent values
    if isempty(opt.hmc_opt) && isempty(opt.ssls_opt) && isempty(opt.sls_opt) && ...
        isempty(opt.latent_opt) && isempty(opt.lik_hmc_opt) && isempty(opt.lik_sls_opt) && ...
        isempty(opt.lik_gibbs_opt)
      opt.ssls_opt.latent_opt.repeat = 20;
      if opt.display>0
        fprintf(' Using SSLS sampler for hyperparameters and ESLS for latent values\n')
      end
    end
    % Set latent values
    if (~isfield(gp,'latentValues') || isempty(gp.latentValues)) ...
        && ~isfield(gp.lik.fh,'trcov')
      if (~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'Softmax', 'Multinom', ...
          'LGP', 'LGPC'}))
        gp.latentValues=zeros(size(y));
      else
        if ~isfield(gp, 'comp_cf') || isempty(gp.comp_cf)
          error('Define multiple covariance functions for latent processes using gp.comp_cf (see gp_set)');
        end        
        if isfield(gp.lik,'xtime')
          ntime = size(gp.lik.xtime,1);
          gp.latentValues=zeros(size(y,1)+ntime,1);
        else
          gp.latentValues=zeros(size(y,1)*length(gp.comp_cf),1);
        end        
      end
    end
  else
    % latent method is not MCMC
    % If no options structures, use SLS as a default sampler for parameters
    if ~isempty(opt.ssls_opt)
      warning('Latent method is not MCMC. ssls_opt ignored')
      opt.ssls_opt=[];
    end
    if ~isempty(opt.latent_opt)
      warning('Latent method is not MCMC. latent_opt ignored')
      opt.latent_opt=[];
    end
    if isempty(opt.hmc_opt) && isempty(opt.sls_opt) && ...
        isempty(opt.lik_hmc_opt) && isempty(opt.lik_sls_opt) && ...
        isempty(opt.lik_gibbs_opt)
      opt.sls_opt.nomit = 0;
      opt.sls_opt.display = 0;
      opt.sls_opt.method = 'minmax';
      opt.sls_opt.wsize = 10;
      opt.sls_opt.plimit = 5;
      opt.sls_opt.unimodal = 0;
      opt.sls_opt.mmlimits = [-10; 10];
      if opt.display>0
        if isfield(gp,'latent_method')
          fprintf(' Using SLS sampler for hyperparameters and %s for latent values\n',gp.latent_method)
        else
          fprintf(' Using SLS sampler for hyperparameters\n')
        end
      end
    end
  end
    
  % Initialize record
  if isempty(record)
    % No old record
    record=recappend();
  else
    ri=size(record.etr,1);
  end

  % Set the states of samplers
  if ~isempty(opt.latent_opt)
    f=gp.latentValues;
    if isfield(opt.latent_opt, 'rstate')
      if ~isempty(opt.latent_opt.rstate)
        latent_rstate = opt.latent_opt.rstate;
      else
        hmc2('state', sum(100*clock))
        latent_rstate=hmc2('state');
      end
    else
      hmc2('state', sum(100*clock))
      latent_rstate=hmc2('state');
    end
  else
    f=y;
  end
  if ~isempty(opt.hmc_opt)
    if isfield(opt.hmc_opt, 'nuts') && opt.hmc_opt.nuts
      % Number of step-size adapting stept in hmc_nuts
      if ~isfield(opt.hmc_opt, 'nadapt')
        opt.hmc_opt.nadapt = 20;
      end
    end
    if isfield(opt.hmc_opt, 'rstate')
      if ~isempty(opt.hmc_opt.rstate)
        hmc_rstate = opt.hmc_opt.rstate;
      else
        hmc2('state', sum(100*clock))
        hmc_rstate=hmc2('state');
      end
    else
      hmc2('state', sum(100*clock))
      hmc_rstate=hmc2('state');
    end
  end    
  if ~isempty(opt.ssls_opt)
    f=gp.latentValues;
  end
  if ~isempty(opt.lik_hmc_opt)
    if isfield(opt.lik_hmc_opt, 'rstate')
      if ~isempty(opt.lik_hmc_opt.rstate)
        lik_hmc_rstate = opt.lik_hmc_opt.rstate;
      else
        hmc2('state', sum(100*clock))
        lik_hmc_rstate=hmc2('state');
      end
    else
      hmc2('state', sum(100*clock))
      lik_hmc_rstate=hmc2('state');
    end        
  end
  
  % Print labels for sampling information
  if opt.display
    fprintf(' cycle  etr      ');
    if ~isempty(opt.hmc_opt)
      fprintf('hrej     ')              % rejection rate of latent value sampling
    end
    if ~isempty(opt.sls_opt)
      fprintf('slsrej  ');
    end
    if ~isempty(opt.lik_hmc_opt)
      fprintf('likel.rej  ');
    end
    if ~isempty(opt.latent_opt)
      if isequal(gp.fh.mc, @esls)
        fprintf('lslsn')              % No rejection rate for esls, print first accepted value
      else
        fprintf('lrej ')              % rejection rate of latent value sampling
      end
      if isfield(opt.latent_opt, 'sample_latent_scale') 
        fprintf('    lvScale    ')
      end
    end
    fprintf('\n');
  end


  % --- Start sampling ------------
  for k=1:opt.nsamples
    
    if opt.persistence_reset
      if ~isempty(opt.hmc_opt)
        hmc_rstate.mom = [];
      end
      if ~isempty(opt.latent_opt)
        if isfield(opt.latent_opt, 'rstate')
          opt.latent_opt.rstate.mom = [];
        end
      end
      if ~isempty(opt.lik_hmc_opt)
        lik_hmc_rstate.mom = [];
      end
    end
    
    hmcrej = 0;
    lik_hmcrej = 0;
    lrej=0;
    indrej=0;
    for l=1:opt.repeat
      
      % --- Sample latent Values  -------------
      if ~isempty(opt.latent_opt)
        [f, energ, diagnl] = gp.fh.mc(f, opt.latent_opt, gp, x, y, z);
        gp.latentValues = f(:);
        f = f(:);
        if ~isequal(gp.fh.mc, @esls)
          lrej=lrej+diagnl.rej/opt.repeat;
        else
          lrej = diagnl.rej;
        end
        if isfield(diagnl, 'opt')
          opt.latent_opt = diagnl.opt;
        end
      end
      
      % --- Sample parameters with HMC ------------- 
      if ~isempty(opt.hmc_opt)
        if isfield(opt.hmc_opt, 'nuts') && opt.hmc_opt.nuts
          % Use NUTS hmc
          w = gp_pak(gp);
          lp = @(w) deal(-gpmc_e(w,gp,x,y,f,z), -gpmc_g(w,gp,x,y,f,z));
          if k<opt.hmc_opt.nadapt
            % Take one sample while adjusting step length
            opt.hmc_opt.Madapt = 1; 
            opt.hmc_opt.M = 0; 
          else
            % Take one sample without adjusting step length
            opt.hmc_opt.Madapt = 0; 
            opt.hmc_opt.M = 1; 
          end
          [w, energies, diagnh] = hmc_nuts(lp, w, opt.hmc_opt);
          opt.hmc_opt = diagnh.opt;
          hmcrej=hmcrej+diagnh.rej/opt.repeat;
          w=w(end,:);
          gp = gp_unpak(gp, w);

        else
          if isfield(opt.hmc_opt,'infer_params')
            infer_params = gp.infer_params;
            gp.infer_params = opt.hmc_opt.infer_params;
          end
          w = gp_pak(gp);
          % Set the state
          hmc2('state',hmc_rstate);
          % sample (y is passed as z, to allow sampling of likelihood parameters)
          [w, energies, diagnh] = hmc2(@gpmc_e, w, opt.hmc_opt, @gpmc_g, gp, x, y, f, z);
          % Save the current state
          hmc_rstate=hmc2('state');
          hmcrej=hmcrej+diagnh.rej/opt.repeat;
          if isfield(diagnh, 'opt')
            opt.hmc_opt = diagnh.opt;
          end
          opt.hmc_opt.rstate = hmc_rstate;
          w=w(end,:);
          gp = gp_unpak(gp, w);
          if isfield(opt.hmc_opt,'infer_params')
            gp.infer_params = infer_params;
          end
        end
      end
      
      % --- Sample parameters with SLS ------------- 
      if ~isempty(opt.sls_opt)
        if isfield(opt.sls_opt,'infer_params')
          infer_params = gp.infer_params;
          gp.infer_params = opt.sls_opt.infer_params;
        end
        w = gp_pak(gp);
        [w, energies, diagns] = sls(@gpmc_e, w, opt.sls_opt, @gpmc_g, gp, x, y, f, z);
        if isfield(diagns, 'opt')
          opt.sls_opt = diagns.opt;
        end
        w=w(end,:);
        gp = gp_unpak(gp, w);
        if isfield(opt.sls_opt,'infer_params')
          gp.infer_params = infer_params;
        end
      end

      % Sample parameters & latent values with SSLS
      if ~isempty(opt.ssls_opt)
        if isfield(opt.ssls_opt,'infer_params')
          infer_params = gp.infer_params;
          gp.infer_params = opt.sls_opt.infer_params;
        end
        w = gp_pak(gp);
        [w, f, diagns] = surrogate_sls(f, w, opt.ssls_opt, gp, x, y, z);
        gp.latentValues = f;
        if isfield(diagns, 'opt')
          opt.ssls_opt = diagns.opt;
        end
        w=w(end,:);
        gp = gp_unpak(gp, w);
        if isfield(opt.sls_opt,'infer_params')
          gp.infer_params = infer_params;
        end
      end      
      
      % --- Sample the likelihood parameters with Gibbs ------------- 
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ...
          isfield(gp.lik,'gibbs') && isequal(gp.lik.gibbs,'on')
        [gp.lik, f] = gp.lik.fh.gibbs(gp, gp.lik, x, f);
      end
      
      % --- Sample the likelihood parameters with HMC ------------- 
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ...
          ~isempty(opt.lik_hmc_opt)
        infer_params = gp.infer_params;
        gp.infer_params = 'likelihood';
        w = gp_pak(gp);
        fe = @(w, lik) (-lik.fh.ll(feval(lik.fh.unpak,lik,w),y,f,z)-lik.fh.lp(feval(lik.fh.unpak,lik,w)));
        fg = @(w, lik) (-lik.fh.llg(feval(lik.fh.unpak,lik,w),y,f,'param',z)-lik.fh.lpg(feval(lik.fh.unpak,lik,w)));
        % Set the state
        hmc2('state',lik_hmc_rstate);
        [w, energies, diagnh] = hmc2(fe, w, opt.lik_hmc_opt, fg, gp.lik);
        % Save the current state
        lik_hmc_rstate=hmc2('state');
        lik_hmcrej=lik_hmcrej+diagnh.rej/opt.repeat;
        if isfield(diagnh, 'opt')
          opt.lik_hmc_opt = diagnh.opt;
        end
        opt.lik_hmc_opt.rstate = lik_hmc_rstate;
        w=w(end,:);
        gp = gp_unpak(gp, w);
        gp.infer_params = infer_params;
      end        
      
      % --- Sample the likelihood parameters with SLS ------------- 
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ...
          ~isempty(opt.lik_sls_opt)
        w = gp_pak(gp, 'likelihood');
        fe = @(w, lik) (-lik.fh.ll(feval(lik.fh.unpak,lik,w),y,f,z) -lik.fh.lp(feval(lik.fh.unpak,lik,w)));
        [w, energies, diagns] = sls(fe, w, opt.lik_sls_opt, [], gp.lik);
        if isfield(diagns, 'opt')
          opt.lik_sls_opt = diagns.opt;
        end
        w=w(end,:);
        gp = gp_unpak(gp, w, 'likelihood');
      end
      
    end % ----- for l=1:opt.repeat ---------  
    
    % --- Set record -------    
    ri=ri+1;
    record=recappend(record);
    
    % Display some statistics  THIS COULD BE DONE NICER ALSO...
    if opt.display && rem(ri,opt.display)==0
      fprintf(' %4d  %.3f  ',ri, record.etr(ri,1));
      if ~isempty(opt.hmc_opt)
        fprintf(' %.1e  ',record.hmcrejects(ri));
      end
      if ~isempty(opt.sls_opt)
        fprintf('sls  ');
      end
      if ~isempty(opt.lik_hmc_opt)
        fprintf(' %.1e  ',record.lik_hmcrejects(ri));
      end
      if ~isempty(opt.latent_opt)
        fprintf('%.1e',record.lrejects(ri));
        fprintf('  ');
        if isfield(diagnl, 'lvs')
          fprintf('%.6f', diagnl.lvs);
        end
      end
      fprintf('\n');
    end
  end
  
%------------------------
function record = recappend(record)
% RECAPPEND - Record append
%          Description
%          RECORD = RECAPPEND(RECORD, RI, GP, P, T, PP, TT, REJS, U) takes
%          old record RECORD, record index RI, training data P, target
%          data T, test data PP, test target TT and rejections
%          REJS. RECAPPEND returns a structure RECORD containing following
%          record fields of:
  
  ncf = length(gp.cf);
  
  if nargin == 0   % Initialize record structure
    record.type = gp.type;
    record.lik = gp.lik;
    if isfield(gp,'latent_method')
      record.latent_method = gp.latent_method;
    end
    if isfield(gp, 'comp_cf')
      record.comp_cf = gp.comp_cf;
    end
    % If sparse model is used save the information about which
    switch gp.type
      case 'FIC'
        record.X_u = [];
      case {'PIC' 'PIC_BLOCK'}
        record.X_u = [];
        record.tr_index = gp.tr_index;
      case 'CS+FIC'
        record.X_u = [];
      otherwise
        % Do nothing
    end
    if isfield(gp,'latentValues')
      record.latentValues = [];
      record.lrejects = 0;
    end
    record.jitterSigma2 = [];
    
    if isfield(gp, 'site_tau')
      record.site_tau = [];
      record.site_nu = [];
      record.Ef = [];
      record.Varf = [];
      record.p1 = [];
    end
    
    % Initialize the records of covariance functions
    for i=1:ncf
      cf = gp.cf{i};
      record.cf{i} = cf.fh.recappend([], gp.cf{i});
      % Initialize metric structure
      if isfield(cf,'metric')
        record.cf{i}.metric = cf.metric.fh.recappend(cf.metric, 1);
      end
    end
    
    % Initialize the record for likelihood
    lik = gp.lik;
    record.lik = lik.fh.recappend([], gp.lik);
    
    % Set the meanfunctions into record if they exist
    if isfield(gp, 'meanf')
      record.meanf = gp.meanf; 
    end
    
    if isfield(gp, 'comp_cf')
      record.comp_cf = gp.comp_cf; 
    end

    if isfield(gp,'p')
      record.p = gp.p;
    end
    if isfield(gp,'latent_method')
      record.latent_method = gp.latent_method;
    end
    if isfield(gp,'latent_opt')
      record.latent_opt = gp.latent_opt;
    end
    if isfield(gp,'fh')
      record.fh=gp.fh;
    end
    
    record.infer_params = gp.infer_params;
    record.e = [];
    record.edata = [];
    record.eprior = [];
    record.etr = [];
    record.hmcrejects = 0;
    ri = 1;
    lrej = 0;
    indrej = 0;
    hmcrej=0;
    lik_hmcrej=0;
  end

  % Set the record for every covariance function
  for i=1:ncf
    gpcf = gp.cf{i};
    record.cf{i} = gpcf.fh.recappend(record.cf{i}, ri, gpcf);
    % Record metric structure
    if isfield(gpcf,'metric')
      record.cf{i}.metric = record.cf{i}.metric.fh.recappend(record.cf{i}.metric, ri, gpcf.metric);
    end
  end

  % Set the record for likelihood
  lik = gp.lik;
  record.lik = lik.fh.recappend(record.lik, ri, lik);

  % Set jitterSigma2 to record
  if ~isempty(gp.jitterSigma2)
    record.jitterSigma2(ri,:) = gp.jitterSigma2;
  end

  % Set the latent values to record structure
  if isfield(gp, 'latentValues')
    record.latentValues(ri,:)=gp.latentValues(:)';
  end

  % Set the inducing inputs in the record structure
  switch gp.type
    case {'FIC', 'PIC', 'PIC_BLOCK', 'CS+FIC'}
      record.X_u(ri,:) = gp.X_u(:)';
  end

  % Record training error and rejects
  if isfield(gp,'latentValues')
    elik = gp.lik.fh.ll(gp.lik, y, gp.latentValues, z);
    [record.e(ri,:),record.edata(ri,:),record.eprior(ri,:)] = gp_e(gp_pak(gp), gp, x, gp.latentValues);
    record.etr(ri,:) = record.e(ri,:) - elik;    
    % Set rejects 
    record.lrejects(ri,1)=lrej;
  else
    [record.e(ri,:),record.edata(ri,:),record.eprior(ri,:)] = gp_e(gp_pak(gp), gp, x, y, 'z', z);
    record.etr(ri,:) = record.e(ri,:);
  end
  
  if ~isempty(opt.hmc_opt)
    record.hmcrejects(ri,1)=hmcrej; 
  end

  if ~isempty(opt.lik_hmc_opt)
    record.lik_hmcrejects(ri,1)=lik_hmcrej; 
  end

  % If inputs are sampled set the record which are on at this moment
  if isfield(gp,'inputii')
    record.inputii(ri,:)=gp.inputii;
  end
  
  if isfield(gp, 'meanf')
      nmf = numel(gp.meanf);
      for i=1:nmf
          gpmf = gp.meanf{i};
          record.meanf{i} = gpmf.fh.recappend(record.meanf{i}, ri, gpmf);
      end
  end
end

function e = gpmc_e(w, gp, x, y, f, z)

  e=0;
  if ~isempty(strfind(gp.infer_params, 'covariance'))
    e=e+gp_e(w, gp, x, f, 'z', z);
  end
  if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
      && ~isfield(gp.lik.fh,'trcov') ...
      && isfield(gp.lik.fh,'lp') && ~isequal(y,f)
    % Evaluate the contribution to the error from non-Gaussian likelihood
    % if latent method is MCMC
    gp=gp_unpak(gp,w);
    lik=gp.lik;
    e=e-lik.fh.ll(lik,y,f,z)-lik.fh.lp(lik);
  end
 
end

function g = gpmc_g(w, gp, x, y, f, z)

  g=[];
  if ~isempty(strfind(gp.infer_params, 'covariance'))
    g=[g gp_g(w, gp, x, f, 'z', z)];
  end
  if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
      && ~isfield(gp.lik.fh,'trcov') ...
      && isfield(gp.lik.fh,'lp') && ~isequal(y,f)
    % Evaluate the contribution to the gradient from non-Gaussian likelihood
    % if latent method is not MCMC
    gp=gp_unpak(gp,w);
    lik=gp.lik;
    g=[g -lik.fh.llg(lik,y,f,'param',z)-lik.fh.lpg(lik)];
  end

end
end
