%DEMO_REGRESSION_ROBUST A regression demo with Student-t
%                       distribution as a residual model.
%
%  Description
%    The synthetic data used here is the same used by Radford M. 
%    Neal in his regression problem with outliers example in
%    Software for Flexible Bayesian Modeling
%    (http://www.cs.toronto.edu/~radford/fbm.software.html). The
%    problem consist of one dimensional input and target variables. 
%    The input data, x, is sampled from standard Gaussian
%    distribution and the corresponding target values come from a
%    distribution with mean given by
%
%           f = 0.3 + 0.4x + 0.5sin(2.7x) + 1.1/(1+x^2).
%
%    For most of the cases the distribution about this mean is
%    Gaussian with standard deviation of 0.1, but with probability
%    0.05 a case is an outlier for wchich the standard deviation is
%    1.0. There are total 200 cases from which the first 100 are
%    used for training and the last 100 for testing.
%
%    We use Student-t distribution as an abservation model
%
%          y ~ St(f, nu, s^2),
%
%    where f is the mean, nu the degrees of freedom and s^2 the
%    scale. The mean is given a GP prior
%
%          f ~ N(0, K).
%
%    The model can be inferred with MCMC or Laplace approximation. 
%    The MCMC can be performed either by utilizing the scale
%    mixture representation of the Student-t distribution or the
%    actual distribution. The scale mixture representation is given
%    as in Gelman et.al (2004)
%
%          y_i ~ N (f_i, a^2*U_i)
%          U_i ~ Inv-Chi^2 (nu, t^2),
%
%    where nu represents the degrees of freedom and a*t = s in the
%    Student-t distribution.
%
%    The demo is organized as follows:
%
%     1) Optimization approach with Normal noise
%     2) MCMC approach with scale mixture noise model (~=Student-t)
%        All parameters sampled
%     3) Laplace approximation Student-t likelihood
%        All parameters optimized
%     4) MCMC approach with Student-t likelihood nu kept fixed to 4
%     5) Laplace approximation Student-t likelihood
%        nu kept fixed to 4
%     6) EP approximation with fractional EP Student-t likelihood
%
%  See Vanhatalo et.al. for discussion on the model and methods.
%
%  Refernces:
%    Vanhatalo, J., Jyl�nki P. and Vehtari, A. (2009). Gaussian
%    process regression with Student-t likelihood. Advances in
%    Neural Information Processing systems
%
%    Gelman, Carlin, Stern and Rubin (2004) Bayesian Data Analysis,
%    second edition. Chapman & Hall / CRC.
%

% Copyright (c) 2010 Jarno Vanhatalo, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% ========================================
% Optimization approach with Normal noise
% ========================================

% load the data. First 100 variables are for training
% and last 100 for test
S = which('demo_regression_robust');
L = strrep(S,'demo_regression_robust.m','demodata/odata.txt');
x = load(L);
y = x(1:100,2);
x = x(1:100,1);
[n, nin] = size(x); 

% Test data
xt = [-2.7:0.01:2.7]';
yt = 0.3+0.4*xt+0.5*sin(2.7*xt)+1.1./(1+xt.^2);

% We create a Gaussian process and priors for GP parameters. Prior for GP
% parameters is Gaussian multivariate hierarchical. The residual is given at
% first Gaussian prior to find good starting value for noise variances

% Construct the priors for the parameters of covariance functions
pl = prior_t();
pm = prior_sqrtunif();
pn = prior_logunif();

% create the Gaussian process
gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
lik = lik_gaussian('sigma2', 0.2^2, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9)

% --- MAP estimate ---
disp('Gaussian noise model and MAP estimate for parameters')

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Prediction
[Eft, Varft, lpyt, Eyt, Varyt] = gp_pred(gp, x, y, xt, 'yt', ones(size(xt)));
std_ft = sqrt(Varft);

% Plot the prediction and data
% plot the training data with dots and the underlying 
% mean of it as a line
figure
hold on
h1=plot(xt,yt, 'k');
h2=plot(xt, Eft, xt, Eft-2*std_ft, 'r--', xt, Eft+2*std_ft, 'r--');
h3=plot(x,y,'b.');
%plot(xt,yt,'r.')
legend([h1 h2(1) h2(3) h3],'real f', 'Ef', 'Ef+-2*std(f)','y',4)
axis on;
title('The predictions and the data points (Gaussian noise model with MAP estimate for parameters)');
drawnow
S1 = sprintf('length-scale: %.3f, magnSigma2: %.3f  \n', ...
             gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2)

% ========================================
% MCMC approach with scale mixture noise model (~=Student-t)
% Here we sample all the variables 
%     (lenghtScale, magnSigma, sigma(noise-t) and nu)
% ========================================
disp(['Scale mixture Gaussian (~=Student-t) noise model                ';...
      'using MCMC integration over the latent values and parameters    '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% Here, set own Sigma2 for every data point
lik = lik_gaussiansmt('ndata', n, 'sigma2', repmat(0.2^2,n,1), ...
                      'nu_prior', prior_logunif());
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);

% Sample
[r,g,opt]=gp_mc(gp, x, y, 'nsamples', 300, 'display', 20);

% thin the record
rr = thin(r,100,2);

figure 
subplot(2,2,1)
hist(rr.lik.nu,20)
title('Mixture model, \nu')
subplot(2,2,2)
hist(sqrt(rr.lik.tau2).*rr.lik.alpha,20)
title('Mixture model, \sigma')
subplot(2,2,3) 
hist(rr.cf{1}.lengthScale,20)
title('Mixture model, length-scale')
subplot(2,2,4) 
hist(rr.cf{1}.magnSigma2,20)
title('Mixture model, magnSigma2')

% make predictions for test set
[Eft, Varft] = gp_pred(rr,x,y,xt);
std_ft=sqrt(Varft);

% Plot the network outputs as '.', and underlying mean with '--'
figure
hold on
h1=plot(xt,yt, 'k');
h2=plot(xt, Eft, xt, Eft-2*std_ft, 'r--', xt, Eft+2*std_ft, 'r--');
h3=plot(x,y,'b.');
%plot(xt,yt,'r.')
legend([h1 h2(1) h2(3) h3],'real f', 'Ef', 'Ef+-2*std(f)','y',4)
axis on;
title('The predictions and the data points (Scale mixture noise model with MCMC)')
drawnow
S2 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', ...
             mean(rr.cf{1}.lengthScale), mean(rr.cf{1}.magnSigma2))

% ========================================
% Laplace approximation Student-t likelihood
%  Here we optimize all the variables 
%  (lengthScale, magnSigma2, sigma(noise-t) and nu)
% ========================================
disp(['Student-t noise model using Laplace integration over the '; ...
      'latent values and MAP estimate for the parameters        '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
pn = prior_logunif();
lik = lik_t('nu', 4, 'nu_prior', pn, ...
            'sigma2', 0.2^2, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-6, ...
            'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Predictions to test points
[Eft, Varft] = gp_pred(gp, x, y, xt);
std_ft = sqrt(Varft);

% Plot the prediction and data
figure
hold on
h1=plot(xt,yt, 'k');
h2=plot(xt, Eft, xt, Eft-2*std_ft, 'r--', xt, Eft+2*std_ft, 'r--');
h3=plot(x,y,'b.');
%plot(xt,yt,'r.')
legend([h1 h2(1) h2(3) h3],'real f', 'Ef', 'Ef+-2*std(f)','y',4)
axis on;
title(sprintf('The predictions and the data points (Student-t noise model (nu=%.2f,sigma2=%.3f) noise) with Laplace+MAP, ',gp.lik.nu, gp.lik.sigma2));
drawnow
S3 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', ...
             gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2)

 
% ========================================
% EP approximation Student-t likelihood
%  Here we optimize all the variables 
%   (lengthScale, magnSigma2, sigma(noise-t) and nu)
% ========================================
disp(['Student-t noise model using EP integration over the';...
      'latent values and MAP estimate for parameters      '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% Create the likelihood structure
pnu = prior_loglogunif();
lik = lik_t('nu', 4, 'nu_prior', pnu, 'sigma2', 0.2^2, ...
            'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4, ...
            'latent_method', 'EP');

% --- MAP estimate ---

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Predictions to test points
[Eft, Varft] = gp_pred(gp, x, y, xt);
std_ft = sqrt(Varft);

% Plot the prediction and data
figure
plot(xt,yt,'k')
hold on
plot(xt,Eft)
plot(xt, Eft-2*std_ft, 'r--')
plot(x,y,'.')
legend('real f', 'Ef', 'Ef+-2*std(f)','y',4)
plot(xt, Eft+2*std_ft, 'r--')
title(sprintf('The predictions and the data points (Student-t noise model, (nu=%.2f,sigma=%.3f) with EP+MAP)',gp.lik.nu, sqrt(gp.lik.sigma2)));
drawnow
S4 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2)
           
           
           
% ========================================
% MCMC approach with Student-t likelihood
%  Here we analyze the model with fixed degrees of freedom
%   nu = 4 
% ========================================
disp(['Student-t noise model with nu= 4 and using MCMC integration';...
      'over the latent values and parameters                      '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
pn = prior_logunif();
lik = lik_t('nu', 4, 'nu_prior', [], 'sigma2', 0.2^2, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4, ...
             'latent_method', 'MCMC');
f=gp_pred(gp,x,y,x);
gp=gp_set(gp, 'latent_opt', struct('f', f));

% Sample 
[rgp,g,opt]=gp_mc(gp, x, y, 'nsamples', 400, 'display', 20);
rr = thin(rgp,100,2);

% make predictions for test set
[Eft, Varft] = gp_pred(rr,x,y,xt);
std_ft = sqrt(Varft);

% Plot the network outputs as '.', and underlying mean with '--'
figure
plot(xt,yt,'k')
hold on
plot(xt,Eft)
plot(xt, Eft-2*std_ft, 'r--')
plot(x,y,'.')
legend('real f', 'Ef', 'Ef+-2*std(f)','y',4)
plot(xt, Eft+2*std_ft, 'r--')
title('The predictions and the data points (Student-t noise model, nu fixed (nu=4), with MCMC)')
drawnow
S5 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', mean(rr.cf{1}.lengthScale), mean(rr.cf{1}.magnSigma2))

% ========================================
% Laplace approximation Student-t likelihood
%  Here we analyze the model with fixed degrees of freedom
%   nu = 4 
% ========================================
disp(['Student-t noise model with nu=4 using Laplace integration over';...
      'the latent values and MAP estimate for the parameters         '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% Create the likelihood structure
pn = prior_logunif();
lik = lik_t('nu', 4, 'nu_prior', [], 'sigma2', 0.2^2, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9, ...
            'latent_method', 'Laplace');
e = gp_e(gp_pak(gp), gp, x, y);

% --- MAP estimate ---

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Predictions to test points
[Eft, Varft] = gp_pred(gp, x, y, xt);
std_ft = sqrt(Varft);

% Plot the prediction and data
figure
plot(xt,yt,'k')
hold on
plot(xt,Eft)
plot(xt, Eft-2*std_ft, 'r--')
plot(x,y,'.')
legend('real f', 'Ef', 'Ef+-2*std(f)','y',4)
plot(xt, Eft+2*std_ft, 'r--')
title(sprintf('The predictions and the data points (Student-t noise model, nu fixed (nu=%.2f,sigma=%.3f) with Laplace+MAP)',gp.lik.nu, sqrt(gp.lik.sigma2)));
drawnow
S6 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2)


% ========================================
% EP approximation Student-t likelihood
%  Here we analyze the model with fixed degrees of freedom
%   nu = 4 
% ========================================
disp(['Student-t noise model with nu=4 using EP integration over';...
      'the latent values and MAP estimate for parameters        '])

gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% Create the likelihood structure
lik = lik_t('nu', 4, 'nu_prior', [], 'sigma2', 0.2^2, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4, ...
            'latent_method', 'EP');

% --- MAP estimate ---

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Predictions to test points
[Eft, Varft] = gp_pred(gp, x, y, xt);
std_ft = sqrt(Varft);

% Plot the prediction and data
figure
plot(xt,yt,'k')
hold on
plot(xt,Eft)
plot(xt, Eft-2*std_ft, 'r--')
plot(x,y,'.')
legend('real f', 'Ef', 'Ef+-2*std(f)','y',4)
plot(xt, Eft+2*std_ft, 'r--')
title(sprintf('The predictions and the data points (Student-t noise model, nu fixed (nu=%.2f,sigma=%.3f) with EP+MAP)',gp.lik.nu, sqrt(gp.lik.sigma2)));
drawnow
S7 = sprintf('length-scale: %.3f, magnSigma2: %.3f \n', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2)
