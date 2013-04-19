%DEMO_PERIODIC  Regression problem demonstration for periodic data
%                    
%  Description
%    A demonstration of the use of periodic covariance function
%    gpcf_periodic with two data sets, the Mauna Loa CO2 data (see
%    for example Rasmussen and Williams 2006) and the monthly
%    Finnish drowning statistics 2002-2008.
%    
%    For the Mauna Loa data, the standard Gaussian process
%    regression model is constructed as in DEMO_REGRESSION2 with
%    additive covariance functions and Gaussian noise
%
%    The observations y are assumed to satisfy
%
%         y = f + g + e,    where e ~ N(0, s^2).
%
%    f and g are underlying latent functions, which we are
%    interested in. We place a zero mean Gaussian process prior
%    them, which implies that at the observed input locations
%    latent values have prior
%
%         f ~ N(0, Kf) and g ~ N(0,Kg)
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters.
%
%    Since both the likelihood and prior are Gaussian, we obtain a
%    Gaussian marginal likelihood
%
%        p(y|th) = N(0, Kf + Kg + I*s^2).
%    
%    By placing a prior for parameters, p(th), we can
%    find the maximum a posterior (MAP) estimate for them by
%    maximizing
%
%       argmax   log p(y|th) + log p(th).
%         th
%   
%    If we want to find an approximation for the posterior of the
%    parameters, we can sample them using Markov chain Monte
%    Carlo (MCMC) methods.
%
%    After finding MAP estimate or posterior samples of
%    parameters, we can use them to make predictions for f.
%
%    For more detailed discussion of Gaussian process regression
%    see, for example, Rasmussen and Williams (2006).
%
%    For the drowning data, a different approach is needed as the
%    likelihood is no longer Gaussian. The regression of counts is
%    implemented with a Poisson likelihood model with Expectation
%    Propagation as the latent optimisation method.
%
%    For details on the implementation see GP_E, GP_G, GP_PRED for
%    the standard regression and GPEP_E, GPEP_G AND GPEP_PRED for
%    the expectation propagation.
%
%  See also  DEMO_REGRESSION2, DEMO_SPATIAL1
%
%  References:
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.

% Copyright (c) 2009 Heikki Peura
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


% This file is organised in two parts:
%  1) Mauna Loa data analysis with GP regression
%  2) Drowning data analysis with Poisson likelihood

%========================================================
% PART 1 Mauna Loa data analysis with full GP model
%========================================================

% Load the data
S = which('demo_periodic');
L = strrep(S,'demo_periodic.m','demodata/maunaloa_data.txt');

data=load(L);
y = data(:, 2:13);
y=y';
y=y(:);
x = [1:1:length(y)]';
x = x(y>0);
y = y(y>0);
avgy = mean(y);
y = y-avgy;

[n,nin] = size(x);
% Now 'x' consist of the inputs and 'y' of the output. 
% 'n' and 'nin' are the number of data points and the 
% dimensionality of 'x' (the number of inputs).

% First, we will do the inference without the periodic covariance function
% (as in DEMO_REGRESSION2), then add the periodic term and compare the
% results

% ---------------------------
% --- Construct the model ---
% 
% First create squared exponential covariance function with ARD and 
% Gaussian noise structures...
gpcf1 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 3);
gpcf2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1);
lik = lik_gaussian();

% ... Then set the prior for the parameters of covariance functions...
pl = prior_t('s2', 3);
pm = prior_sqrtunif();
gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp(gpcf2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', {gpcf1,gpcf2});

% -----------------------------
% --- Conduct the inference ---
%
% We will make the inference first by finding a maximum a posterior estimate 
% for the parameters via gradient based optimization.  

% --- MAP estimate ---
% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-4);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Make predictions. Below Eyt_full is the predictive mean and
% Varyt_full the predictive variance.
xt=[1:650]';

[Eft_full, Varft_full, lpyt_full, Eyt_full, Varyt_full] = gp_pred(gp, x, y, xt, 'yt', ones(650,1));

% Plot the prediction and data
figure;hold on
plot(x,y,'.', 'MarkerSize',7)
plot(xt,Eyt_full,'k', 'LineWidth', 2)
plot(xt,Eyt_full-2.*sqrt(Varyt_full),'k--')
plot(xt,Eyt_full+2.*sqrt(Varyt_full),'k--')
axis tight
caption1 = sprintf('GP with sexp+sexp+noise:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f \n s^2_{noise} = %.2f', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2, gp.cf{2}.lengthScale, gp.cf{2}.magnSigma2, gp.lik.sigma2);
title(caption1)
legend('Data point', 'predicted mean', '2\sigma error', 'Location', 'NorthWest')

% -------------------------------------------
% INFERENCE WITH PERIODIC COVARIANCE FUNCTION

% With the increasing number of parameters, the optimisation takes
% longer, especially with period length optimisation included. The results
% are however significantly better. Both models fit the data well, yet only
% the one with the periodic component has real predictive power.

% ---------------------------
% --- Construct the model ---
% 
% First create a set of covariance functions: a long term squared
% exponential function, two short term ones, the periodic function and a
% noise structure
gpcf1 = gpcf_sexp('lengthScale', 67*12, 'magnSigma2', 66*66);
gpcfp = gpcf_periodic('lengthScale', 1.3, 'magnSigma2', 2.4*2.4);
gpcfp = gpcf_periodic(gpcfp, 'period', 12,'lengthScale_sexp', 90*12, 'decay', 1);
lik = lik_gaussian('sigma2', 0.3);
gpcf2 = gpcf_sexp('lengthScale', 2, 'magnSigma2', 2);

% ... Then set the prior for the parameters of covariance functions...
pl = prior_t('s2', 10, 'nu', 3);
pn = prior_t('s2', 10, 'nu', 4);

gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pl);
gpcf2 = gpcf_sexp(gpcf2, 'lengthScale_prior', pl, 'magnSigma2_prior', pl);
gpcfp = gpcf_periodic(gpcfp, 'lengthScale_prior', pl, 'magnSigma2_prior', pl);
gpcfp = gpcf_periodic(gpcfp, 'lengthScale_sexp_prior', pl, 'period_prior', pn);
lik = lik_gaussian(lik, 'sigma2_prior', pn);

% ... Finally create the GP structure
gp = gp_set('lik', lik, 'cf', {gpcf1, gpcfp, gpcf2});

% -----------------------------
% --- Conduct the inference ---
%
% We will make the inference first by finding a maximum a posterior
% estimate for the parameters via gradient based optimization.

% --- MAP estimate ---
% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-4,'Display','iter');
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Make predictions. Below Eft_full is the predictive mean and
% Varft_full the predictive variance.
xt=[1:650]';

[Eft_full, Varft_full, lpyt_full, Eyt_full, Varyt_full] = gp_pred(gp, x, y, xt, 'yt', ones(650,1));

% Plot the prediction and data
figure;hold on
plot(xt,Eyt_full,'k')
plot(xt,Eyt_full-2.*sqrt(Varyt_full),'k--')
plot(xt,Eyt_full+2.*sqrt(Varyt_full),'k--')
plot(x,y,'.', 'MarkerSize',7)
axis tight
caption1 = sprintf('GP sexp+periodic+sexp+noise:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f, p=%.2f, s_sexp^2 = %.2f, \n l_3= %.2f, s^2_3 = %.2f, \n l_4= %.2f, s^2_4 = %.2f, \n s^2_{noise} = %.2f', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2, gp.cf{2}.lengthScale, gp.cf{2}.magnSigma2, gp.cf{2}.period, gp.cf{2}.lengthScale_sexp, gp.cf{3}.lengthScale, gp.cf{3}.magnSigma2, gp.lik.sigma2);
title(caption1)
legend('Data point', 'predicted mean', '2\sigma error','Location','NorthWest')

% Plot the latent components separately
[Eft_full1, Varft_full1] = gp_pred(gp, x, y, x, 'predcf', 1);
[Eft_full2, Varft_full2] = gp_pred(gp, x, y, x, 'predcf', [2 3]);

figure
[AX, H1, H2] = plotyy(x, Eft_full2, x, Eft_full1);
set(H2,'LineStyle','--')
set(H2, 'LineWidth', 2)
%set(H1, 'Color', 'k')
set(H1,'LineStyle','-')
set(H1, 'LineWidth', 0.8)
title('The long and short term latent component')

%========================================================
% PART 2 Drowning data analysis with FULL GP
%========================================================

% Here we use a GP model with Poisson likelihood to analyse the
% monthly Finnish drowning mortality data from 2002-2008. Finland,
% with almost 200 000 lakes and a long coast on the Baltic sea, has a
% relatively high drowning mortality among the developed countries.
% It is well known that drownings exhibit a periodic behaviour within
% the year, peaking in the summer holiday season in July and coming to
% near zero in the winter when most lakes and the Baltic sea are
% frozen.

% The Poisson likelihood is chosen to deal with the regression of
% counts.  As the amount of drownings, although small in the
% wintertime, can never be negative, a Gaussian likelihood is not
% suitable. A negative binomial is another option, especially with
% overdispersed data (ie. with high variance), as it provides another
% parameter to control the dispersion.

% Load the data

S = which('demo_periodic');
L = strrep(S,'demo_periodic.m','demodata/drowning.txt');
data=load(L);
y = data(:, 2:13);
y=y';
y=y(:);
y1=y;
y=y(1:72);
avgy = mean(y);
x = [1:length(y)]';

[n,nin] = size(x);

% ---------------------------
% --- Construct the model ---
% 

% Create covariance functions. Here we use a squared exponential
% and a neural network function to deal with long term change,
% another SE for short term effects and a periodic component for
% the cyclic nature of the data. The period of the cycle is not
% optimised as it is strongly believed to be exactly 12 months.

gpcf1 = gpcf_sexp('lengthScale', [67], 'magnSigma2', 1);
gpcfp = gpcf_periodic('lengthScale', [1.3], 'magnSigma2', 2.4*2.4,...
    'period', 12,'lengthScale_sexp', 50, 'decay', 1);
likn=gpcf_neuralnetwork('biasSigma2',10,'weightSigma2',3);
gpcf2 = gpcf_sexp('lengthScale', [2], 'magnSigma2', 2);

% ... Then set the prior for the parameters of covariance functions...
pl = prior_t('s2', 1000, 'nu', 3);
pm = prior_sqrtt('s2', 2, 'nu', 3);
pl2 = prior_t('s2', 5, 'nu', 3);
pm2 = prior_sqrtt('s2', 3, 'nu', 3);
ppl = prior_t('s2', 100, 'nu', 3);
ppm = prior_sqrtt('s2', 1, 'nu', 3);
pn = prior_t('s2', 10, 'nu', 4);
ppp = prior_t('s2', 100, 'nu', 4);

gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp(gpcf2, 'lengthScale_prior', pl2, 'magnSigma2_prior', pm2);
likn = gpcf_neuralnetwork(likn, 'biasSigma2_prior', pn, 'weightSigma2_prior', ppp);
gpcfp = gpcf_periodic(gpcfp, 'lengthScale_prior', ppl, 'magnSigma2_prior', ppm,  'lengthScale_sexp_prior', pl);
lik = lik_gaussian(lik, 'sigma2_prior', pn);

% ... Create the GP structure, Poisson likelihood with
% Expectation Propagation as approximation method
z=repmat(mean(y),length(y),1);
gp = gp_set('lik', lik_poisson(), 'cf', {gpcf1,gpcfp,gpcf2,likn});
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-4);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'z', z,'opt',opt,'optimf',@fminlbfgs);

% Predictions
xt=[1:96]';
xt=[1:132]';
[Eft_full, Varft_full] = gp_pred(gp, x, y, xt, 'z', z);

% Plot results
xtt=2001+23/24+1/12*xt;
figure;hold on
plot(xtt(1:length(y),1),y,'b.', 'MarkerSize',20)
plot(xtt(length(y)+1:length(y1),1),y1(length(y)+1:length(y1)),'b*', 'MarkerSize',7)
plot(xtt(:,1),exp(Eft_full).*mean(y),'b', 'LineWidth', 2)
plot(xtt(:,1),exp(Eft_full-1.96.*sqrt(Varft_full)).*mean(y),'b--')
plot(xtt(:,1),exp(Eft_full+1.96.*sqrt(Varft_full)).*mean(y),'b--')

legend('Training data', 'Validation data','Predicted mean','95% CI', 'Location', 'NorthWest')
line(2008,0:80,'LineWidth',2)
axis tight
