%DEMO_MODELASSESMENT2  Demonstration for model assessment when the observation 
%                      model is non-Gaussian
%
%  Description
%    We will consider the classification problem in demo_classific. 
%    The analysis is conducted with full Gaussian process using
%    both probit and logit likelihood. The performance of these two
%    models are compared by evaluating the ten-fold
%    cross-validation, leave-one-out cross-validation, WAIC, DIC
%    and the effective number of parameters The inference will be
%    conducted using maximum a posterior (MAP) estimate for the
%    parameters using EP and Laplace approximation, via full Markov
%    chain Monte Carlo (MCMC) and with an integration approximation
%    (IA) for the parameters.
%
%    This demo is organised in two parts:
%     1) data analysis with with probit likelihood
%     2) data analysis with with logit likelihood
%
%  See also  
%    DEMO_CLASSIFIC1, DEMO_MODELASSESMENT1

% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


% =====================================
% 1) data analysis with probit likelihood
% =====================================
disp('Data analysis with probit likelihood')
S = which('demo_classific');
L = strrep(S,'demo_classific.m','demodata/synth.tr');
x=load(L);
y=x(:,end);
y = 2.*y-1;
x(:,end)=[];
[n, nin] = size(x);

DIC=repmat(NaN,1,8);DIC2=repmat(NaN,1,8);DIC_latent=repmat(NaN,1,8);
p_eff=repmat(NaN,1,8);p_eff2=repmat(NaN,1,8);p_eff_latent=repmat(NaN,1,8);p_eff_latent2=repmat(NaN,1,8);

% Create covariance functions
gpcf = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 2);

% Set the prior for the parameters of covariance functions 
pl = prior_logunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl,'magnSigma2_prior', pl); %

% Create the GP structure
gp = gp_set('lik', lik_probit, 'cf', gpcf, 'jitterSigma2', 1e-4);

% ------- Laplace approximation --------
disp(['Probit with Laplace integration over the latent values '; ...
      'and MAP estimate for the parameters                    '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'Laplace');

n=length(y);

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{1} = 'pr_Laplace';
p_eff_latent(1) = gp_peff(gp, x, y);
[DIC_latent(1), p_eff_latent2(1)] = gp_dic(gp, x, y, 'focus', 'latent');
WAIC(1) = gp_waic(gp,x,y);

% Evaluate the 10-fold cross-validation results. 
cvres = gp_kfcv(gp, x, y, 'display', 'fold');
mlpd_cv(1) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp, x, y);
mlpd_loo(1) = mean(lpy);

% ------- Expectation propagation --------
disp(['Probit with EP integration over the latent values and MAP '; ...
      'estimate for the parameters                               '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{2} = 'pr_EP';
p_eff_latent(2) = gp_peff(gp, x, y) ;
[DIC_latent(2), p_eff_latent2(2)] = gp_dic(gp, x, y, 'focus', 'latent');
WAIC(2) = gp_waic(gp,x,y);

% Evaluate the 10-fold cross-validation results. 
cvres = gp_kfcv(gp, x, y, 'display', 'fold');
mlpd_cv(2) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp, x, y);
mlpd_loo(2) = mean(lpy);

% ------- MCMC ---------------
disp(['Probit with MCMC integration over the latent values and '; ...
      'the parameters                                          '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'MCMC');

% Sample
mcopt.nsamples=220;mcopt.display=20;
[rgp,gp]=gp_mc(gp, x, y, mcopt);
rgp=thin(rgp, 21, 2);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{3} = 'pr_MCMC';
[DIC(3), p_eff(3)] =  gp_dic(rgp, x, y, 'focus', 'param');
[DIC2(3), p_eff2(3)] =  gp_dic(rgp, x, y, 'focus', 'all');
WAIC(3) = gp_waic(rgp,x,y);

% Evaluate the 10-fold cross-validation results. 
mcopt.nsamples=50;mcopt.display=20;
cvres = gp_kfcv(gp, x, y, 'inf_method', 'MCMC', 'opt', mcopt, 'display', 'fold');
mlpd_cv(3) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(rgp, x, y);
mlpd_loo(3) = mean(lpy);

% --- Integration approximation approach ---
disp(['Probit with EP integration over the latent values and '; ...
      'grid integration over the parameters                  '])

% Use EP
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% now perform the integration
clear opt
opt.int_method = 'grid';
opt.step_size = 2;
gp_array = gp_ia(gp, x, y, opt);

models{4} = 'pr_IA'; 
[DIC(4), p_eff(4)] =  gp_dic(gp_array, x, y, 'focus', 'param');
[DIC2(4), p_eff2(4)] =  gp_dic(gp_array, x, y, 'focus', 'all');
WAIC(4) = gp_waic(gp_array,x,y);

% Then the 10 fold cross-validation.
cvres = gp_kfcv(gp, x, y, 'inf_method', 'IA', 'opt', opt, 'display', 'fold');
mlpd_cv(4) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp_array, x, y);
mlpd_loo(4) = mean(lpy);

% =====================================
% 2) data analysis with logit likelihood
% =====================================
disp('Data analysis with logit likelihood')

S = which('demo_classific');
L = strrep(S,'demo_classific.m','demodata/synth.tr');
x=load(L);
y=x(:,end);
y = 2.*y-1;
x(:,end)=[];
[n, nin] = size(x);

% Create covariance functions
gpcf = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 2);

% Set the prior for the parameters of covariance functions 
pl = prior_logunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl,'magnSigma2_prior', pl); %

% Create the likelihood structure
lik = ('init');

% Create the GP structure
gp = gp_set('lik', lik_logit, 'cf', gpcf, 'jitterSigma2', 1e-4);


% ------- Laplace approximation --------
disp(['Logit with Laplace integration over the latent values and '; ...
      'MAP estimate for the parameters                           '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{5} = 'lo_Laplace';
p_eff_latent(5) = gp_peff(gp, x, y);
[DIC_latent(5), p_eff_latent2(5)] = gp_dic(gp, x, y, 'focus', 'latent');
WAIC(5) = gp_waic(gp,x,y);

% Evaluate the 10-fold cross-validation results. 
cvres = gp_kfcv(gp, x, y, 'display', 'fold');
mlpd_cv(5) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp, x, y);
mlpd_loo(5) = mean(lpy);

% ------- Expectation propagation --------
disp(['Logit with EP integration over the latent values and MAP'; ...
      'estimate for the parameters                             '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{6} = 'lo_EP';
p_eff_latent(6) = gp_peff(gp, x, y) ;
[DIC_latent(6), p_eff_latent2(6)] = gp_dic(gp, x, y, 'focus', 'latent');
WAIC(6) = gp_waic(gp,x,y);

% Evaluate the 10-fold cross-validation results. 
cvres = gp_kfcv(gp, x, y, 'display', 'fold');
mlpd_cv(6) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp, x, y);
mlpd_loo(6) = mean(lpy);

% ------- MCMC ---------------
disp(['Logit with MCMC integration over the latent values and '; ...
      'the parameters                                         '])

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'MCMC');

% Sample 
mcopt.nsamples=200;mcopt.display=20;
[rgp,gp] = gp_mc(gp, x, y, mcopt);
rgp=thin(rgp, 21, 2);

% Evaluate the effective number of parameters and DIC with focus on latent variables.
models{7} = 'lo_MCMC';
[DIC(7), p_eff(7)] =  gp_dic(rgp, x, y, 'focus', 'param');
[DIC2(7), p_eff2(7)] =  gp_dic(rgp, x, y, 'focus', 'all');
WAIC(7) = gp_waic(rgp,x,y);

% Evaluate the 10-fold cross-validation results. 
mcopt.nsamples=50;mcopt.display=20;
cvres = gp_kfcv(gp, x, y, 'inf_method', 'MCMC', 'opt', mcopt, 'display', 'fold');
mlpd_cv(7) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(rgp, x, y);
mlpd_loo(7) = mean(lpy);

% --- Integration approximation approach ---
disp(['Logit with EP integration over the latent values and grid '; ...
      'integration over the parameters                           '])

% Use EP
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% now perform the integration
clear opt
opt.int_method = 'grid';
opt.step_size = 2;
gp_array = gp_ia(gp, x, y, opt);

models{8} = 'lo_IA'; 
[DIC(8), p_eff(8)] =  gp_dic(gp_array, x, y, 'focus', 'param');
[DIC2(8), p_eff2(8)] =  gp_dic(gp_array, x, y, 'focus', 'all');
WAIC(8) = gp_waic(gp_array,x,y);

% Then the 10 fold cross-validation.
cvres = gp_kfcv(gp, x, y, 'inf_method', 'IA', 'opt', opt, 'display', 'fold');
mlpd_cv(8) = cvres.mlpd_cv;

% Evaluate the leave-one-out cross-validation results. 
[Ef,Varf,lpy] =  gp_loopred(gp_array, x, y);
mlpd_loo(8) = mean(lpy);

%========================================================
% PART 4 Print the results
%========================================================
disp('Summary of the results')

S = '       ';
for i = 1:length(models)
    S = [S '  ' models{i}];
end

S = sprintf([S '\n CV-mlpd  %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], mlpd_cv);
S = sprintf([S '\n LOO-mlpd %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], mlpd_loo);
S = sprintf([S '\n ']);
S = sprintf([S '\n WAIC     %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], WAIC);
S = sprintf([S '\n ']);
S = sprintf([S '\n DIC_h    %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], DIC);
S = sprintf([S '\n DIC_a    %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], DIC2);
S = sprintf([S '\n DIC_l    %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], DIC_latent);
S = sprintf([S '\n peff_h   %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], p_eff);
S = sprintf([S '\n peff_a   %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], p_eff2);
S = sprintf([S '\n peff_l   %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], p_eff_latent);
S = sprintf([S '\n peff_l2  %6.2f    %6.2f  %6.2f  %6.2f   %6.2f    %6.2f  %6.2f  %6.2f'], p_eff_latent2);
S = sprintf([S '\n ']);
S = sprintf([S '\n The notation is as follows:']);
S = sprintf([S '\n pr_*     = probit likelihood and inference method']);
S = sprintf([S '\n lo_*     = logit likelihood and inference method']);
S = sprintf([S '\n CV-mlpd  = mean log predictive density from the 10-fold CV. ']);
S = sprintf([S '\n LOO-mlpd = mean log predictive density from the 10-fold CV. ']);
S = sprintf([S '\n WAIC     = Widely applicable information criterion. ']);
S = sprintf([S '\n DIC_h    = DIC with focus on parameters. ']);
S = sprintf([S '\n DIC_a    = DIC with focus on parameters and laten variables (all). ']);
S = sprintf([S '\n DIC_l    = DIC with focus on latent variables. ']);
S = sprintf([S '\n peff_h   = effective number of parameters (latent variables marginalized). ']);
S = sprintf([S '\n peff_a   = effective number of parameters and latent variables. ']);
S = sprintf([S '\n peff_l   = effective number of latent variables evaluated with gp_peff. ']);
S = sprintf([S '\n peff_l2  = effective number of latent variables evaluated with gp_dic. ']);
S = sprintf([S '\n '])
