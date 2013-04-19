function [criteria, cvpreds, cvws, trpreds, trw, cvtrpreds] = gp_kfcv(gp, x, y, varargin)
%GP_KFCV  K-fold cross validation for GP model
%
%  Description
%    [CRITERIA, CVPREDS, CVWS, TRPREDS, TRW] = GP_KFCV(GP, X, Y, OPTIONS)
%    Performs K-fold cross-validation for a GP model given input matrix X
%    and target vector Y.
%
%    OPTIONS is optional parameter-value pair
%      z          - optional observed quantity in triplet (x_i,y_i,z_i)
%                   Some likelihoods may use this. For example, in
%                   case of Poisson likelihood we have z_i=E_i,
%                   that is, expected value for ith case.
%      yt         - optional observed yt in test points. 
%                   Default option for yt is yt=y. Note: if given yt
%                   criteria will be empty.
%      pred       - optional return values for cvpreds (f,lp,y). 
%                   Default option for pred is 'f+lp+y'
%      inf_method - inference method. Possible methods are
%                    'MAP'      parameters optimized to MAP (default)
%                    'MCMC'     MCMC sampling using GP_MC
%                    'IA'       integration approximation using GP_IA
%                    'fixed'    parameters are fixed, it either use MAP 
%                               or integration approximation, depending if 
%                               GP is a single GP structure or a GP array
%                               (for example from GP_IA)
%      optimf     - function handle for an optimization function, which is
%                   assumed to have similar input and output arguments
%                   as usual fmin*-functions. Default is @fminscg.
%      opt        - options for the inference method. If 'MAP' is used
%                   use optimset to set options for optimization.
%                   Default options for optimization are 'GradObj'
%                   is 'on', 'LargeScale' is 'off', 'Display' is 'off'
%      k          - number of folds in CV, default k=10
%      rstream    - number of a random stream to be used for
%                   permuting the data befor division. This way
%                   same permutation can be obtained for different
%                   models. Default is 1. See doc RandStream for
%                   more information.
%      trindex    - k-fold CV training indices. A cell array with k
%                   fields each containing index vector for respective
%                   training set.
%      tstindex   - k-fold CV test indices. A cell array with k
%                   fields each containing index vector for
%                   respective test set.
%      display    - defines if messages are displayed. 
%                   - 'off' displays no output
%                   - 'on' (default) gives some output  
%                   - 'iter' displays output at each iteration
%      save_results
%                 - defines if detailed results are stored 'false'
%                   (default) or 'true'. If 'true' gp_kfcv stores the
%                   results in the current working directory into a
%                   cv_resultsX folder (or in 'folder', see next
%                   option), where X is a number. If there are
%                   cv_results* folders already, X is the smallest
%                   number not in use yet.
%       folder    - string defining the folder where to save the
%                   results. That is, the results will be stored in
%                   'current working directory'/folder. See previous
%                   option for default.
%
%    The output arguments are the following
%       criteria  - structure including the following fields
%                     mlpd_cv     - mean log predictive density
%                     Var_lpd_cv  - variance estimate for mlpd
%                     rmse_cv     - root mean squared error
%                     Var_rmse_cv - variance estimate for mrmse
%                     mabs_cv     - mean absolute error
%                     Var_abs_cv  - variance estimate for mabs
%       cvpreds   - CV predictions structure including the same fields
%                   as trpreds
%       trpreds   - training predictions structure including
%                   the following fields:
%                     Eft
%                     Varft
%                     Eyt
%                     Varyt
%                     pyt
%       cvws      - parameter weight vectors for each CV fold
%       trw       - parameter weight vector for training data
%       cvtrpreds - training and test predictions structure for 
%                   every cv fold including the following fields:
%                     Eft
%                     Varft
%                     Eyt
%                     Varyt
%                     pyt
%
%     The K-fold cross validation is performed as follows: The data
%     are divided into k groups D_k. For each group, we evaluate
%     the test statistics
%
%            u(D_k | D_{k-1})
%
%     where u is the utility function and D_{k-1} is the data in
%     the k-1 groups other than k. The utility functions provided
%     by gp_kfcv are
%
%       log predictive density
%         lpd(D_k | D_{k-1})  = mean( log( p( y_k|D_{k-1} ) ) )
%       squared error
%         rmse(D_k | D_{k-1}) = mean( ( E[y_k|D_{k-1}] - y_k ).^2 )
%       absolute error
%         abs(D_k | D_{k-1})  = mean( abs( E[y_k|D_{k-1}] - y_k ) )
%
%     After the utility is evaluated for each group, we can
%     evaluate the output arguments, which are obtained as follows
%
%       mean log predictive density
%         mlpd_cv  = mean( lpd(D_k | D_{k-1}) )          ,k=1...K
%       root mean squared error
%         mrmse_cv = sqrt( mean( rmse(D_k | D_{k-1}) ) ) ,k=1...K
%       mean absolute error
%         mabs_cv  = mean( abs(D_k | D_{k-1}) )          ,k=1...K
%
%     The variance estimates for the above statistics are evaluated
%     across the groups K. For mean log predictive density and mean
%     absolute error this reduces to evaluate, for example,
%
%         Var_lpd_cv = var( lpd(D_k | D_{k-1}) ) / K,    k=1...K.
%
%     For root mean squared error, we need to take the square root
%     of each group statistics first to obtain
%
%         Var_rmse_cv = var( sqrt( rmse(D_k | D_{k-1}) ) ) / K,    k=1...K.
%
%     The above statistics are returned by the function. However,
%     if we use the save_results option we obtain some additional
%     test statistics, which are only saved in the result file.
%     These extra statistics include, for example, bias corrected
%     expected utilities (Vehtari and Lampinen, 2002) and the
%     training utility for the whole data and each cross-validation
%     training set. The detailed list of variables saved in the
%     result file is:
%
%     For more information see the file cv_results.mat,
%     which contains the following variables
%       lpd_cv      - log predictive density (nx1 vector)
%       rmse_cv     - squared error (nx1 vector)
%       abs_cv      - absolute error (nx1 vector)
%       mlpd_cv     - mean log predictive density (a scalar summary)
%       mrmse_cv    - root mean squared error (a scalar summary)
%       mabs_cv     - mean absolute error (a scalar summary)
%       Var_lpd_cv  - variance of mean log predictive density
%                     (a scalar summary)
%       Var_rmse_cv - variance of the root mean squared error
%                     (a scalar summary)
%       Var_abs_cv  - variance of the mean absolute error
%                     (a scalar summary)
%       trindex     - training indices
%       tstindex    - test indices
%       lpd_cvtr    - mean log predictive density for each of
%                     k-CV training sets (kx1 vector)
%       rmse_cvtr   - root mean squared error for each of
%                     k-CV training sets (kx1 vector)
%       abs_cvtr    - absolute error for each of
%                     k-CV training sets (kx1 vector)
%       lpd_tr      - log predictive density for the
%                     full training set
%       rmse_tr     - root mean squared error for the
%                     full training set
%       abs_tr      - absolute error for the full trainng set
%       mlpd_ccv    - log predictive density with corrected k-CV
%       rmse_ccv    - root mean squared error with corrected k-CV
%       abs_ccv     - absolute error with corrected k-CV
%       cpu_time    - the cpu time used for inferring the full
%                     data set
%
%  See also
%    DEMO_MODELASSESMENT1, GP_PEFF, GP_DIC
%
%  References:
%    Spiegelhalter, Best, Carlin and van der Linde (2002). Bayesian
%    measures of model complexity and fit. J. R. Statist. Soc. B,
%    64, 583-639.
%
%    Gelman, Carlin, Stern and Rubin (2004) Bayesian Data Analysis,
%    second edition. Chapman & Hall / CRC.
%
%    Aki Vehtari and Jouko Lampinen. Bayesian model assessment and
%    comparison using cross-validation predictive densities. Neural
%    Computation, 14(10):2439-2468, 2002.
%

%  Experimental features
%      inf_method - inference method. Possible methods are
%                    'LOO'      parameters optimized using leave-one-out
%                    'KFCV'     parameters optimized using k-fold-CV
%                    'WAIC'     parameters optimized using WAIC
  
% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari
% Copyright (c) 2012 Ernesto Ulloa

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_KFCV';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('inf_method', 'MAP', @(x) ...
                   ismember(x,{'MAP' 'LOO' 'KFCV' 'WAIC' 'WAICV' 'WAICG' 'MCMC' 'IA' 'fixed'}))
  ip.addParamValue('optimf', @fminscg, @(x) isa(x,'function_handle'))
  ip.addParamValue('opt', struct(), @isstruct)
  ip.addParamValue('pred','f+lp+y',@ischar)
  ip.addParamValue('k', 10, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
  ip.addParamValue('rstream', 1, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
  ip.addParamValue('trindex', [], @(x) isempty(x) || iscell(x))
  ip.addParamValue('tstindex', [], @(x) isempty(x) || iscell(x))
  ip.addParamValue('display', 'on', @(x) islogical(x) || ...
                   ismember(x,{'on' 'off' 'iter' 'fold'}))
  ip.addParamValue('save_results', false, @(x) islogical(x))
  ip.addParamValue('folder', [], @(x) ischar(x) )
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;
  yt=ip.Results.yt;
  inf_method=ip.Results.inf_method;
  optimf=ip.Results.optimf;
  opt=ip.Results.opt;
  pred=ip.Results.pred;
  k=ip.Results.k;
  rstream=ip.Results.rstream;
  trindex=ip.Results.trindex;
  tstindex=ip.Results.tstindex;
  display = ip.Results.display;
  if isequal(display,'fold');display='iter';end
  save_results=ip.Results.save_results;
  folder = ip.Results.folder;

  [n,nin] = size(x);
  gp_orig = gp;

  if ismember(inf_method,{'MAP' 'LOO' 'KFCV' 'WAIC' 'WAICV' 'WAICG'})
    optdefault=struct('Display','off');
    opt=optimset(optdefault,opt);
  end

  if (isempty(trindex) && ~isempty(tstindex)) || (~isempty(trindex) && isempty(tstindex))
    error('gp_kfcv: If you give cross-validation indexes, you need to provide both trindex and tstindex.')
  end

  if isempty(trindex) || isempty(tstindex)
    [trindex, tstindex] = cvit(n, k, rstream);
  else
    k=length(trindex);
  end
  parent_folder = pwd;

  cvws=[];
  trw=[];
  % loop over the crossvalidation sets
  if ismember(display,{'iter'})
    fprintf('Evaluating the CV utility. The inference method is %s.\n',inf_method)
  end
  
  % *** note: yt must be a scalar or a vector of size y
  
  if ~isempty(yt)
    if size(yt,1)==size(y,1)
      yt=yt;
    elseif size(yt,1)==1
      yt=bsxfun(@times,ones(size(y)),yt);
    else
      error('size of yt does not match y nor is it a scalar');   
    end
    ytflag=1;
  else
    ytflag=0;
    yt=y;
  end
  
  if ~isempty(strfind(pred, 'f'))
    predft=1;
  else
    predft=0;
  end
  
  if ~isempty(strfind(pred, 'lp'))
    predlpyt=1;
  else 
    predlpyt=0;
  end
  
  if ~isempty(strfind(pred, 'y'))
    predyt=1;
  else
    predyt=0;
  end
  
  % Initialize empty structures
  Eft=[]; Varft=[]; lpyt=[]; Eyt=[]; Varyt=[];
  cvpreds.Eft=[]; cvpreds.Varft=[]; cvpreds.lpyt=[];
  cvpreds.Eyt=[]; cvpreds.Varyt=[];
  cvtrpreds.Eft=[]; cvtrpreds.Varft=[]; cvtrpreds.lpyt=[];
  cvtrpreds.Eyt=[]; cvtrpreds.Varyt=[];
  lpd_cv =[]; lpd_cvtr = []; lpd_cvm = [];
  rmse_cv = []; rmse_cvtr = []; abs_cv= []; abs_cvtr = [];  
  rmse_cvm = []; abs_cvm = [];
  mlpd_cv = []; Var_lpd_cv = []; mrmse_cv = [];
  mabs_cv = []; Var_rmse_cv = []; Var_abs_cv = [];
  criteria.mlpd_cv=[]; criteria.Var_lpd_cv=[];
  criteria.mrmse_cv=[]; criteria.Var_rmse_cv=[];
  criteria.mabs_cv=[]; criteria.Var_abs_cv=[];
  
  
  for i=1:k
    if isempty(tstindex{i})
      continue
    end

    if isequal(display,'iter')
      fprintf(' The CV-fold number: %d/%d \n', i, k)
    end

    % Set the training and test sets for i'th cross-validation set
    xtr = x(trindex{i},:);
    ytr = y(trindex{i},:);
    yttr= yt(trindex{i},:);
    xtst = x(tstindex{i},:);
    ytst = y(tstindex{i},:);
    yttst = yt(tstindex{i},:);

    if ~isempty(z)
      ztr = z(trindex{i},:);
      zt = z;
      %* yt = y;
      %       opt_tst.zt = z;
      %       opt_tst.yt = y;
    else
      ztr = [];
      %*yt = y;
      zt = [];
      %       opt_tr = struct();
      %       opt_tst.yt = y;
    end

    gp = gp_orig;

    if iscell(gp)
      gptype=gp{1}.type;
    else
      gptype=gp.type;
    end
    tstind=[];
    switch gptype
      case {'FIC' 'CS+FIC'}
        tstind = trindex{i};
      case 'PIC'
        % Set the block indices for the cv set of data points. 
        ntr = size(xtr,1);
        ntst = size(xtst,1);
        trind = [];
        for i1=1:length(gp.tr_index)
          tstind{i1} = [];
          trind{i1} = [];
          for j1 = 1:length(gp.tr_index{i1})
            indtmp = find(all(bsxfun(@minus,xtr,x(gp.tr_index{i1}(j1),:))==0,2));
            %find( sum((xtr - repmat(x(gp.tr_index{i1}(j1),:),ntr,1)).^2,2) == 0 );
            if isempty( indtmp )
              indtmp = find(all(bsxfun(@minus,xtst,x(gp.tr_index{i1}(j1),:))==0,2));
              %find( sum((xtst - repmat(x(gp.tr_index{i1}(j1),:),ntst,1)).^2,2) == 0 );
              tstind{i1} = [tstind{i1} indtmp];
            else
              trind{i1} = [trind{i1} indtmp];
            end
          end
        end
        if iscell(gp)
          for j=1:numel(gp)
            gp{j}.tr_index=trind;
          end
        else
          gp.tr_index = trind;
        end
        tstind=gp_orig.tr_index;
    end

    % Conduct inference
    switch inf_method
      case 'MAP'
        gp=gp_optim(gp,xtr,ytr,'z',ztr,'opt',opt, 'optimf', optimf);
        w=gp_pak(gp);
        cvws(i,:)=w;
      case {'LOO' 'KFCV' 'WAIC' 'WAICV' 'WAICG'}
        if ismember('optimf',ip.UsingDefaults)
          gp=gp_optim(gp,xtr,ytr,'z',ztr,'opt',opt,'loss',inf_method);
        else
          gp=gp_optim(gp,xtr,ytr,'z',ztr,'opt',opt,'loss',inf_method, 'optimf', optimf);
        end
        w=gp_pak(gp);
        cvws(i,:)=w;
      case 'MCMC'
        if numel(gp.jitterSigma2)>1
          gp=thin(gp,numel(gp.jitterSigma2)-1);
        end
        % Scaled mixture noise model is a special case
        % where we need to modify the noiseSigmas2 vector
        % to a right length
        if isequal(gp.lik.type, 'lik_smt')
          gp.lik.noiseSigmas2 = gp_orig.lik.noiseSigmas2(trindex{i});
          gp.lik.r = gp_orig.lik.r(trindex{i});
          gp.lik.U = gp_orig.lik.U(trindex{i});
          gp.lik.ndata = length(trindex{i});
        end
        % Pick latent values for the training set in this fold
        if isfield(gp,'latentValues')
          if (~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'Softmax', 'Multinom', ...
              'LGP', 'LGPC'}))
            latentValues=reshape(gp_orig.latentValues, size(y,1), size(y,2));
            gp.latentValues=reshape(latentValues(trindex{i},:), size(y,2)*length(trindex{i}), 1);
            % gp.latentValues=gp_orig.latentValues(trindex{i});
          else
            if ~isfield(gp.lik, 'xtime')
              nl=length(gp.comp_cf);
              gp.latentValues=gp_orig.latentValues(trindex{i}+(0:nl-1)*n);
            else
              ntime=size(gp.lik.xtime,1);
              gp.latentValues=gp_orig.latentValues([1:ntime, (ntime+trindex{i})]);
            end
          end
        end
        gp = gp_mc(gp, xtr, ytr, 'z', ztr, opt);
        nburnin = floor(length(gp.etr)/3);
        gp = thin(gp,nburnin);
      case 'IA'
        gp = gp_ia(gp, xtr, ytr, 'z', ztr, opt);
      case 'fixed'
        % nothing to do here
    end

    if iscell(gp)
      gplik=gp{1}.lik;
    else
      gplik=gp.lik;
    end

    if predyt
      [Eft, Varft, lpyt, Eyt, Varyt] = gp_pred(gp, xtr, ytr, x, 'tstind', ...
                                     tstind, 'z', ztr, 'yt', yt, 'zt', zt);
    elseif predlpyt
      [Eft, Varft,lpyt] = gp_pred(gp, xtr, ytr, x, 'tstind', ...
                        tstind, 'z', ztr, 'yt', yt, 'zt', zt);
    elseif predft
      [Eft,Varft] = gp_pred(gp, xtr, ytr, x, 'tstind', ...
                    tstind, 'z', ztr, 'yt', yt, 'zt', zt);
    end
    if ~predft
      Eft=[]; Varft=[];
    end
    if ~predlpyt
      lpyt=[];
    end


    if predyt && (isempty(Eyt)||isempty(Varyt))
      warning('This likelihood does not return Eyt and/or Varyt. Empty values will be given');
      emptyyt=1; 
    elseif isempty(Eyt)||isempty(Varyt)
      emptyyt=1;
    else
      emptyyt=0;
    end 

    if predft && (isempty(Eft)||isempty(Varft))
      warning('This likelihood does not return Eft and/or Varft. Empty values will be given');
      emptyft=1;
    elseif isempty(Eft)||isempty(Varft)
      emptyft=1;
    else
      emptyft=0;
    end

    if predlpyt && isempty(lpyt)
      warning('This likelihood does not return lpyt. Empty values will be given');
      emptylp=1;
    elseif isempty(lpyt)
      emptylp=1;
    else
      emptylp=0;
    end
    
    % Save results to structures
    
    if nargout>=6

      if predft && ~emptyft
        cvtrpreds.Eft([trindex{i}(:) ; tstindex{i}(:)],i)=Eft([trindex{i}(:) ; tstindex{i}(:)],:);
        cvtrpreds.Varft([trindex{i}(:) ; tstindex{i}(:)],i)=Varft([trindex{i}(:) ; tstindex{i}(:)],:);
      end

      if predlpyt && ~emptylp
        cvtrpreds.lpyt([trindex{i}(:) ; tstindex{i}(:)],i)=lpyt([trindex{i}(:) ; tstindex{i}(:)],:); 
      end

      if predyt && ~emptyyt
        cvtrpreds.Eyt([trindex{i}(:) ; tstindex{i}(:)],i)=Eyt([trindex{i}(:) ; tstindex{i}(:)],:);
        cvtrpreds.Varyt([trindex{i}(:) ; tstindex{i}(:)],i)=Varyt([trindex{i}(:) ; tstindex{i}(:)],:);   
      end

    end
    
    if nargout>=2

      if predft  && ~emptyft
        cvpreds.Eft(tstindex{i},:)=Eft(tstindex{i},:);
        cvpreds.Varft(tstindex{i},:)=Varft(tstindex{i},:);
      end

      if predlpyt && ~emptylp
        cvpreds.lpyt(tstindex{i},:)=lpyt(tstindex{i},:);               
      end

      if predyt && ~emptyyt 
        cvpreds.Eyt(tstindex{i},:)=Eyt(tstindex{i},:);
        cvpreds.Varyt(tstindex{i},:)=Varyt(tstindex{i}(:),:);
      end

    end

    if predlpyt && ~emptylp
      lpd_cv(tstindex{i}) = log(mean(exp(lpyt(tstindex{i},:)),2));
      lpd_cvtr(i) = mean(log(mean(exp(lpyt(trindex{i})),2)));
      lpd_cvm(i) = mean(log(mean(exp(lpyt(tstindex{i},:)),2)));
    end

    if predyt && ~emptyyt 
      rmse_cv(tstindex{i}) = (mean(Eyt(tstindex{i},:),2) - ytst).^2;
      rmse_cvtr(i) = sqrt(mean((mean(Eyt(trindex{i},:),2) - ytr).^2));

      abs_cv(tstindex{i}) = abs(mean(Eyt(tstindex{i},:),2) - ytst);
      abs_cvtr(i) = mean(abs(mean(Eyt(trindex{i},:),2) - ytr));

      rmse_cvm(i) = sqrt(mean((mean(Eyt(tstindex{i},:),2) - ytst).^2));
      abs_cvm(i) = mean(abs(mean(Eyt(tstindex{i},:),2) - ytst)); 

    end

  end

  if predlpyt && ~emptylp
    mlpd_cv = mean(lpd_cv);
    Var_lpd_cv = var(lpd_cvm)./k;
  end
  
  
  if predyt && ~emptyyt
    mrmse_cv = sqrt(mean(rmse_cv));
    mabs_cv = mean(abs_cv);
    Var_rmse_cv = var(rmse_cvm)./k;
    Var_abs_cv = var(abs_cvm)./k; 
  end
  
  % *** note: if user sends yt then criteria should be empty
  if ~ytflag
    if predlpyt
      criteria.mlpd_cv=mlpd_cv;
      criteria.Var_lpd_cv=Var_lpd_cv;
    end
    if ~emptyyt && predyt
      criteria.mrmse_cv=mrmse_cv;
      criteria.Var_rmse_cv=Var_rmse_cv;
      criteria.mabs_cv=mabs_cv;
      criteria.Var_abs_cv=Var_abs_cv; 
    end 
  end
  
  if save_results || nargout >=4
    % compute full training result

    gp = gp_orig;
    if iscell(gp)
      gptype=gp{1}.type;
    else
      gptype=gp.type;
    end
    switch gptype
      case {'FULL' 'VAR' 'DTC' 'SOR'}
        tstind = [];
      case 'FIC'
        tstind = 1:n;
      case 'PIC'
        if iscell(gp)
          tstind = gp{1}.tr_index;
        else
          tstind = gp.tr_index;
        end
    end
    
    if ~isempty(z)
      opt_tr.z = z;
      opt_tst.zt = z;
      %opt_tst.yt=y;
      opt_tst.yt = yt;
    else
      opt_tr = struct();
      %opt_tst.yt=y;
      opt_tst.yt = yt;
    end
    
    % Evaluate the training utility
    if ismember(display,{'iter'})
      fprintf('\n Evaluating the training utility \n')
    end
    
    % Initialize empty structures
    

    % Conduct inference
    cpu_time = cputime;
    switch inf_method
      case 'MAP'
        gp=gp_optim(gp,x,y,'z',z,'opt',opt_tr, 'optimf', optimf);
        w=gp_pak(gp);
        trw=w;
      case {'LOO' 'KFCV' 'WAIC' 'WAICV' 'WAICG'}
        gp=gp_optim(gp,x,y,'z',z,'opt',opt_tr,'loss',inf_method, 'optimf', optimf);
        w=gp_pak(gp);
        trw=w;
      case 'MCMC'
        gp = gp_mc(gp, x, y, opt_tr, opt);
        nburnin = floor(length(gp.etr)/3);
        gp = thin(gp,nburnin);
      case 'IA'
        gp = gp_ia(gp, x, y, opt_tr, opt);
      case 'fixed'
        % nothing to do here
    end
    cpu_time = cputime - cpu_time;

    % Initialize empty structures
    lpd_tr=[]; rmse_tr=[]; abs_tr=[]; mlpd_ccv=[]; 
    mrmse_ccv=[]; mabs_ccv=[];
    criteria.mlpd_ccv=[]; criteria.rmse_cvv=[]; criteria.mabs_ccv=[];
        
    % Make the predictions
    if predyt
      [Eft, Varft, lpyt, Eyt, Varyt] = gp_pred(gp, x, y, x, 'tstind', tstind, opt_tr, opt_tst); 
    elseif predlpyt
      [Eft, Varft,lpyt] = gp_pred(gp, x, y, x, 'tstind', tstind, opt_tr, opt_tst);
    elseif predft
      [Eft,Varft] = gp_pred(gp, x, y, x, 'tstind', tstind, opt_tr, opt_tst);
    end
    if ~predft
      Eft=[]; Varft=[];
    end
    if ~predlpyt
      lpyt=[];
    end
    
    if predyt && (isempty(Eyt)||isempty(Varyt))
      warning('This likelihood does not return Eyt and/or Varyt empty values will be given');
      emptyyt=1; 
    elseif isempty(Eyt)||isempty(Varyt)
      emptyyt=1;
    else
      emptyyt=0;
    end 

    if predft && (isempty(Eft)||isempty(Varft))
      warning('This likelihood does not return Eft and/or Varft empty values will be given');
      emptyft=1;
    elseif isempty(Eft)||isempty(Varft)
      emptyft=1;
    else
      emptyft=0;
    end

    if predlpyt && isempty(lpyt)
      warning('This likelihood does not return lpyt empty values will be given');
      emptylp=1;
    elseif isempty(lpyt)
      emptylp=1;
    else
      emptylp=0;
    end
    
    % Save predictions to structures
    if nargout>=4
      if predft && ~emptyft  
        trpreds.Eft=Eft;
        trpreds.Varft=Varft;
      end
      if predlpyt && ~emptylp
        trpreds.lpyt=lpyt;
      end
      if predyt && ~emptyyt
        trpreds.Eyt=Eyt;
        trpreds.Varyt=Varyt;
      end
    end

    if predlpyt && ~emptylp
      lpd_tr = mean(log(mean(exp(lpyt),2)));
    end
    if predyt && ~emptyyt
      rmse_tr = sqrt(mean((mean(Eyt,2) - y).^2));
      abs_tr = mean(abs(mean(Eyt,2) - y));
    end

    % compute bias corrected results
    if predlpyt && ~emptylp
      mlpd_ccv =  mlpd_cv +  mean(lpd_tr) -  mean(lpd_cvtr);
    end
    if predyt && ~emptyyt
      mrmse_ccv =  mrmse_cv +  mean(rmse_tr) -  mean(rmse_cvtr); 
      mabs_ccv =  mabs_cv +  mean(abs_tr) -  mean(abs_cvtr);
    end
    
    % *** note: if user sends yt then criteria should be empty 
    if ~ytflag
      if predlpyt && ~emptylp
        criteria.mlpd_ccv=mlpd_ccv;
      end

      if predyt && ~emptyyt
        criteria.rmse_ccv=mrmse_ccv;
        criteria.mabs_ccv=mabs_ccv;
      end
    end

  end

  if save_results
    % Save the results
    if isempty(folder)
      succes = 1;
      result = 1;
      while succes == 1
        folder = sprintf('cv_results%d', result);
        if exist(['./' folder])
          result = result + 1;
        else
          succes = 0;
        end
      end
    else
      if exist(folder)
        folder_alk = folder;
        succes = 1;
        result = 1;
        while succes == 1
          folder = sprintf([folder '%d'], result);
          if exist(['./' folder])
            result = result + 1;
          else
            succes = 0;
          end
        end
        warning('The given folder: %s exists already. gp_kfcv saves the results in: %s instead.', folder_alk, folder)
      end
    end
    mkdir(folder);

    save([folder '/cv_results.mat'], 'lpd_cv', 'rmse_cv', 'abs_cv', 'mlpd_cv', 'mrmse_cv',...
         'mabs_cv','Var_lpd_cv', 'Var_rmse_cv', 'Var_abs_cv', 'trindex', 'tstindex', 'lpd_cvtr', 'rmse_cvtr',...
         'abs_cvtr', 'lpd_tr', 'rmse_tr', 'abs_tr', 'mlpd_ccv', 'mrmse_ccv', 'mabs_ccv', 'cpu_time', 'cvpreds');

    if ismember(display,{'on','iter'})
      fprintf('The results have been saved in the folder:\n %s/%s \n', parent_folder, folder);
    end

    f = fopen([folder '/description.txt'],'w');
    fprintf(f,'The cv results were the following: \n\n');
    fprintf(f,'mlpd_cv = %.4f (+/- %.4f) \n', mlpd_cv, Var_lpd_cv);
    fprintf(f,'mrmse_cv = %.4f (+/- %.4f)  \n', mrmse_cv, Var_rmse_cv);
    fprintf(f,'mabs_cv = %.4f (+/- %.4f) \n', mabs_cv, Var_abs_cv);
    fprintf(f,'cpu time = %.4f \n', cpu_time);
    fprintf(f,'\n\n');
    fprintf(f,'For more information see the file cv_results.mat, which contains ');
    fprintf(f,'the following variables\n\n');
    fprintf(f,'lpd_cv      = log predictive density (nx1 vector) \n');
    fprintf(f,'rmse_cv     = squared error (nx1 vector) \n');
    fprintf(f,'abs_cv      = absolute error (nx1 vector) \n');
    fprintf(f,'mlpd_cv     = mean log predictive density (a scalar summary) \n');
    fprintf(f,'mrmse_cv    = root mean squared error (a scalar summary) \n');
    fprintf(f,'mabs_cv     = mean absolute error (a scalar summary) \n');
    fprintf(f,'Var_lpd_cv  = variance of mean log predictive density (a scalar summary) \n');
    fprintf(f,'Var_rmse_cv = variance of the root mean squared error (a scalar summary) \n');
    fprintf(f,'Var_abs_cv  = variance of the mean absolute error (a scalar summary) \n');
    fprintf(f,'trindex     = training indices \n');
    fprintf(f,'tstindex    = test indices \n');
    fprintf(f,'lpd_cvtr    = mean log predictive density for each of k-CV trainng sets (kx1 vector) \n');
    fprintf(f,'rmse_cvtr   = root mean squared error for each of k-CV trainng sets (kx1 vector) \n');
    fprintf(f,'abs_cvtr    = absolute error for each of k-CV training sets (kx1 vector) \n');
    fprintf(f,'lpd_tr      = log predictive density for the full training set  \n');
    fprintf(f,'rmse_tr     = root mean squared error for the full training set  \n');
    fprintf(f,'abs_tr      = absolute error for the full training set  \n');
    fprintf(f,'lpd_ccv     = log predictive density with corrected cross validation   \n');
    fprintf(f,'rmse_ccv    = root mean squared error with corrected cross validation \n');
    fprintf(f,'abs_ccv     = absolute error with corrected cross validation \n');
    fprintf(f,'cpu_time    = The cpu time used for inferring the full data set \n');
    fprintf(f,'cvpreds     = CV predictions structure  \n');
    fclose(f);
  end
end