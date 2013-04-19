%DEMO_REGRESSION1  Regression problem demonstration for 2-input 
%                  function with Gaussian process
%
%  Description
%    The regression problem consist of a data with two input
%    variables and one output variable with Gaussian noise. The
%    model constructed is following:
%
%    The observations y are assumed to satisfy
%
%         y = f + e,    where e ~ N(0, s^2)
%
%    where f is an underlying function, which we are interested in. 
%    We place a zero mean Gaussian process prior for f, which
%    implies that at the observed input locations latent values
%    have prior
%
%         f ~ N(0, K),
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters.
%
%    Since both likelihood and prior are Gaussian, we obtain a
%    Gaussian marginal likelihood
%
%        p(y|th) = N(0, K + I*s^2).
%    
%   By placing a prior for parameters, p(th), we can find
%   the maximum a posterior (MAP) estimate for them by maximizing
%
%       argmax   log p(y|th) + log p(th).
%         th
%   
%   An approximation for the posterior of the parameters, can be
%   found using Markov chain Monte Carlo (MCMC) methods. We can
%   integrate over the parameters also with other integration
%   approximations such as grid integration.
%
%   After finding MAP estimate or posterior samples of
%   parameters, we can use them to make predictions for f_new:
%
%       p(f_new | y, th) = N(m, S),
%
%          m = K_nt*(K + I*s^2)^(-1)*y
%          S = K_new - K_nt*(K + I*s^2)^(-1)*K_tn
%   
%   where K_new is the covariance matrix of new f, and K_nt between
%   new f and training f.
%
%   For more detailed discussion of Gaussian process regression see,
%   for example, Rasmussen and Williams (2006) or Vanhatalo and
%   Vehtari (2008)
%
%   The demo is organised in three parts:
%     1) data analysis with MAP estimate for the parameters
%     2) data analysis with grid integration over the parameters
%     3) data analysis with MCMC integration over the parameters
%
%  See also DEMO_*
%
%  References:
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%    Vanhatalo, J. and Vehtari, A. (2008). Modelling local and global
%    phenomena with sparse Gaussian processes. Proceedings of the 24th
%    Conference on Uncertainty in Artificial Intelligence,

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

%========================================================
% PART 1 data analysis with full GP model
%========================================================
disp('GP with Gaussian noise model')

% Load the data
S = which('demo_regression1');
L = strrep(S,'demo_regression1.m','demodata/dat.1');
data=load(L);
x = [data(:,1) data(:,2)];
y = data(:,3);
[n, nin] = size(x);

% Now 'x' consist of the inputs and 'y' of the output. 
% 'n' and 'nin' are the number of data points and the 
% dimensionality of 'x' (the number of inputs).

% ---------------------------
% --- Construct the model ---
% 
% First create structures for Gaussian likelihood and squared
% exponential covariance function with ARD
lik = lik_gaussian('sigma2', 0.2^2);
gpcf = gpcf_sexp('lengthScale', [1.1 1.2], 'magnSigma2', 0.2^2)

% Set some priors
pn = prior_logunif();
lik = lik_gaussian(lik,'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Following lines do the same since the default type is FULL
%gp = gp_set('type','FULL','lik',lik,'cf',gpcf);
gp = gp_set('lik', lik, 'cf', gpcf);

% Demostrate how to evaluate covariance matrices. 
% K contains the covariance matrix without noise variance 
%  at the diagonal (the prior covariance)
% C contains the covariance matrix with noise variance at 
% the diagonal (the posterior covariance)
example_x = [-1 -1 ; 0 0 ; 1 1];
[K, C] = gp_trcov(gp, example_x)

% What has happend this far is the following
% - we created structures 'gpcf' and 'lik', which describe 
%   the properties of the covariance function and Gaussian likelihood (see
%   gpcf_sexp and lik_gaussian for more details)
% - we created structures that describe the prior of the length-scale 
%   and magnitude of the squared exponential covariance function and
%   the prior of the noise variance. These structures were set into
%   'gpcf' and 'lik' (see prior_* for more details)
% - we created a GP structure 'gp', which has among others 'gpcf' 
%   and 'lik' structures.  (see gp_set for more details)

% -----------------------------
% --- Conduct the inference ---
%
% We will make the inference first by finding a maximum a posterior
% estimate for the parameters via gradient based optimization. 
% After this we will use grid integration and Markov chain Monte
% Carlo sampling to integrate over the parameters.
 

% --- MAP estimate  ---
disp(' MAP estimate for the parameters')
% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% get optimized parameter values for display
[w,s]=gp_pak(gp);
% display exp(w) and labels
disp(s), disp(exp(w))

% For last, make predictions of the underlying function on a dense
% grid and plot it. Below Eft_map is the predictive mean and
% Varf_map the predictive variance.
[xt1,xt2]=meshgrid(-1.8:0.1:1.8,-1.8:0.1:1.8);
xt=[xt1(:) xt2(:)];
[Eft_map, Varft_map] = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure(1)
clf
mesh(xt1, xt2, reshape(Eft_map,37,37));
hold on
plot3(x(:,1), x(:,2), y, '*')
axis on;
title('The predicted underlying function and the data points (MAP solution)');

% --- Grid integration ---
disp(' Grid integration over the parameters')
% Perform the grid integration and make predictions
% We can also get predictions for marginals of f, which are not
% Gaussian due to integration over the hyperparameters
[gp_array, P_TH, th, Eft_ia, Varft_ia, fx_ia, x_ia] = ...
    gp_ia(gp, x, y, xt, 'int_method', 'grid');

% Plot the predictions for two input locations
figure(2)
clf
subplot(2,1,1)
plot(x_ia(100,:), fx_ia(100,:))
title('p(f|D) at input location (-1.6, 0.7)');
subplot(2,1,2)
plot(x_ia(400,:), fx_ia(400,:))
title('p(f|D) at input location (-0.8, 1.1)');

% --- MCMC ---
disp(' MCMC integration over the parameters')
[gp_rec,g,opt] = gp_mc(gp, x, y, 'nsamples', 220,'display',20);

% After sampling we delete the burn-in and thin the sample chain
gp_rec = thin(gp_rec, 21, 2);

% Make the predictions
[Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);

figure(1)
clf
subplot(1,2,1)
mesh(xt1, xt2, reshape(Eft_map,37,37));
hold on
plot3(x(:,1), x(:,2), y, '*')
axis on;
title(['The predicted underlying function ';
       'and the data points (MAP solution)']);
subplot(1,2,2)
mesh(xt1, xt2, reshape(Eft_mc,37,37));
hold on
plot3(x(:,1), x(:,2), y, '*')
axis on;
title(['The predicted underlying function  ';
       'and the data points (MCMC solution)']);
set(gcf,'pos',[93 511 1098 420])

% We can compare the posterior samples of the parameters to the MAP
% estimate that we got from optimization
figure(3)
clf
subplot(1,2,1)
plot(gp_rec.cf{1}.lengthScale)
title('The sample chain of length-scales')
subplot(1,2,2)
plot(gp_rec.cf{1}.magnSigma2)
title('The sample chain of magnitude')
set(gcf,'pos',[93 511 1098 420])

figure(4)
clf
subplot(1,4,1)
hist(gp_rec.cf{1}.lengthScale(:,1))
hold on
plot(gp.cf{1}.lengthScale(1), 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Length-scale 1')
subplot(1,4,2)
hist(gp_rec.cf{1}.lengthScale(:,2))
hold on
plot(gp.cf{1}.lengthScale(2), 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Length-scale 2')
subplot(1,4,3)
hist(gp_rec.cf{1}.magnSigma2)
hold on
plot(gp.cf{1}.magnSigma2, 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('magnitude')
subplot(1,4,4)
hist(gp_rec.lik.sigma2)
hold on
plot(gp.lik.sigma2, 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Noise variance')
legend('MCMC samples', 'MAP estimate')
set(gcf,'pos',[93 511 1098 420])


% Sample from two posterior marginals and plot them alongside 
% with the MAP and grid integration results
% gpmc_preds returns the predictive mean of the latent function
% with every sampled parameter value.
[Eft_mcs, Varft_mcs] = gpmc_preds(gp_rec, x, y, xt);
sf = normrnd(Eft_mcs(100,:), sqrt(Varft_mcs(100,:)));
sf2 = normrnd(Eft_mcs(400,:), sqrt(Varft_mcs(400,:)));

figure(2)
subplot(1,2,1)
[N,X] = hist(sf);
hist(sf)
hold on
plot(x_ia(100,:), max(N)/max(fx_ia(100,:))*fx_ia(100,:), 'k')
ff = norm_pdf(x_ia(100,:)', Eft_map(100), sqrt(Varft_map(100)));
plot(x_ia(100,:), max(N)/max(ff)*ff, 'r', 'lineWidth', 2)
set(gca, 'Ytick', [])
title('p(f|D) at input location (-1.6, 0.7)');
xlim([0 1])

subplot(1,2,2)
[N,X] = hist(sf2);
hist(sf2)
hold on
plot(x_ia(400,:), max(N)/max(fx_ia(400,:))*fx_ia(400,:), 'k')
ff = norm_pdf(x_ia(400,:)', Eft_map(400), sqrt(Varft_map(400)));
plot(x_ia(400,:), max(N)/max(ff)*ff, 'r', 'lineWidth', 2)
set(gca, 'Ytick', [])
title('p(f|D) at input location (-0.8, 1.1)');
xlim([-1.2 -0.6])

disp('Done')


% ========================
% Print figures for manual
% ========================
% $$$ sf = normrnd(Eft_mc(100,:), sqrt(Varft_mc(100,:)));
% $$$ sf2 = normrnd(Eft_mc(400,:), sqrt(Varft_mc(400,:)));
% $$$ 
% $$$ figure
% $$$ subplot(1,2,1)
% $$$ [N,X] = hist(sf);
% $$$ hist(sf)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(x_ia(100,:), max(N)/max(fx_ia(100,:))*fx_ia(100,:), 'k')
% $$$ ff = norm_pdf(x_ia(100,:)', Eft_map(100), sqrt(Varft_map(100)));
% $$$ plot(x_ia(100,:), max(N)/max(ff)*ff, 'k', 'lineWidth', 2)
% $$$ set(gca, 'Ytick', [])
% $$$ xlim([0 1])
% $$$ ylim([0 110])
% $$$ 
% $$$ subplot(1,2,2)
% $$$ [N,X] = hist(sf2);
% $$$ hist(sf2)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(x_ia(400,:), max(N)/max(fx_ia(400,:))*fx_ia(400,:), 'k')
% $$$ ff = norm_pdf(x_ia(400,:)', Eft_map(400), sqrt(Varft_map(400)));
% $$$ plot(x_ia(400,:), max(N)/max(ff)*ff, 'k', 'lineWidth', 2)
% $$$ set(gca, 'Ytick', [])
% $$$ xlim([-1.2 -0.5])
% $$$ 
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 7 5])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_regression1_fig3.eps
% $$$ 
% $$$ 
% $$$ figure(4)
% $$$ clf, subplot(1,4,1)
% $$$ hist(gp_rec.cf{1}.lengthScale(:,1))
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(gp.cf{1}.lengthScale(1), 0, 'kx', 'MarkerSize', 11, 'LineWidth', 2)
% $$$ xlabel('Length-s 1')
% $$$ xlim([0.3 1.6])
% $$$ 
% $$$ subplot(1,4,2)
% $$$ hist(gp_rec.cf{1}.lengthScale(:,2))
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(gp.cf{1}.lengthScale(2), 0, 'kx', 'MarkerSize', 11, 'LineWidth', 2)
% $$$ xlabel('Length-s 2')
% $$$ xlim([0.4 1.4])
% $$$ 
% $$$ subplot(1,4,3)
% $$$ hist(gp_rec.cf{1}.magnSigma2)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(gp.cf{1}.magnSigma2, 0, 'kx', 'MarkerSize', 11, 'LineWidth', 2)
% $$$ xlabel('magnitude')
% $$$ xlim([0.5 6])
% $$$ 
% $$$ subplot(1,4,4)
% $$$ hist(gp_rec.lik.sigma2)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ plot(gp.lik.sigma2, 0, 'kx', 'MarkerSize', 11, 'LineWidth', 2)
% $$$ xlabel('Noise variance')
% $$$ xlim([0.03 0.06])
% $$$ set(gca, 'Xtick', [0.03 0.06])
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 11 5])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_regression1_fig2.eps
