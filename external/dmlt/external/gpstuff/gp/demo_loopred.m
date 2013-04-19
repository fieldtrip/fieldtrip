%DEMO_LOOPRED  Leave-one-out prediction demonstration for 2 classes
%
%  Description
%    The demonstration program is based on synthetic two class data
%    used by B.D. Ripley (Pattern Recognition and Neural Networks,
%    1996}. The data consists of 2-dimensional vectors that are
%    divided into two classes, labeled 0 or 1. Each class has a
%    bimodal distribution generated from equal mixtures of Gaussian
%    distributions with identical covariance matrices. A Bayesian
%    approach is used to find the decision line and predict the
%    classes of new data points.
%
%    The probability of y being one is assumed to be 
%
%      p(y=1|f) = 1 / (1+exp(-f))
%
%    The latent values f are given a zero mean Gaussian process
%    prior. This implies that at the observed input locations
%    latent values have prior
%
%      f ~ N(0, K),
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters.
% 
%    Here we demonstrate Laplace leave-one-out with 3 different methods:
%    Linear response style, EP-style cavity and Inla-style, EP-LOO and
%    compare LOO-approximations to brute force LOO-CV.
%
%   References:
%
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%    Neal, R. M. (1998) Regression and classification using
%    Gaussian process priors (with discussion), in J. M. Bernardo,
%    et al (editors) Bayesian Statistics 6, Oxford University
%    Press, pp. 475-501:
%

% Copyright (c) 2013 Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% This demonstration is based on the dataset used in the book Pattern
% Recognition and Neural Networks by B.D. Ripley (1996), Cambridge
% University Press.

% Training data
S = which('demo_classific');
L = strrep(S,'demo_classific.m','demodata/synth.tr');
x=load(L);
y=x(:,end);
y = 2.*y-1;
x(:,end)=[];
[n, nin] = size(x);

% Test data
xt1=repmat(linspace(min(x(:,1)),max(x(:,1)),20)',1,20);
xt2=repmat(linspace(min(x(:,2)),max(x(:,2)),20)',1,20)';
xt=[xt1(:) xt2(:)];

% Create likelihood function
lik = lik_probit();
%lik = lik_logit();

% Create covariance functions
gpcf = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 10);

% Set the prior for the parameters of covariance functions 
pl = prior_t();
pm = prior_sqrtunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl,'magnSigma2_prior', pm); %

% Create the GP structure (type is by default FULL)
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);

% ------- Laplace approximation --------
fprintf(['%s model with Laplace integration over the latent values and\n' ...
         'MAP estimate for the parameters\n'],gp.lik.type)

% Set the approximate inference method 
% (Laplace is the default, so this could be skipped)
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Linear-response style leave-one-out (default)
[Eft_lrs, Varft_lrs, lpyt_lrs] = gpla_loopred(gp, x, y, 'method', 'lrs');

% EP-style cavity leave-one-out 
[Eft_cav, Varft_cav, lpyt_cav] = gpla_loopred(gp, x, y, 'method', 'cavity');

% Inla-style leave-one-out 
[Eft_inla, Varft_inla, lpyt_inla] = gpla_loopred(gp, x, y, 'method', 'inla');

% LOO-CV
% Here we should use 250 folds, but for the sake of computation time, only
% 100 folds are used (results are almost the same)
[crit, pred] = gp_kfcv(gp, x, y, 'opt', opt, 'k', 100, 'Display', 'iter');

% ------- EP --------
fprintf(['%s model with EP integration over the latent values and\n' ...
         'MAP estimate for the parameters\n'],gp.lik.type)

% Set the approximate inference method 
% (Laplace is the default, so this could be skipped)
gp = gp_set(gp, 'latent_method', 'EP');

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% EP-LOO
[Eft_ep, Varft_ep, lpyt_ep] = gpep_loopred(gp, x, y);

% LOO-CV
[crit_ep, pred_ep] = gp_kfcv(gp, x, y, 'opt', opt, 'k', 100, 'Display', 'iter');

% Plot results
figure, subplot(2,2,1)
plot(pred.lpyt, lpyt_lrs,'.');ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(lrs) vs. LOO-CV');
xlabel('Log-predictive density (CV)'); ylabel('Log-predictive density(LA-LOO)');

subplot(2,2,2), plot(pred.lpyt, lpyt_cav,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(cavity) vs. LOO-CV');
xlabel('Log-predictive density (CV)'); ylabel('Log-predictive density(LA-LOO)');

subplot(2,2,3), plot(pred.lpyt, lpyt_inla,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(inla) vs. LOO-CV');
xlabel('Log-predictive density (CV)'); ylabel('Log-predictive density(LA-LOO)');

subplot(2,2,4), plot(pred_ep.lpyt, lpyt_ep,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('EP-LOO vs. LOO-CV');
xlabel('Log-predictive density (CV)'); ylabel('Log-predictive density(EP-LOO)');

figure, subplot(2,2,1)
plot(pred.Eft, Eft_lrs,'.');ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(lrs) vs. LOO-CV');
xlabel('Latent prediction (CV)'); ylabel('Latent prediction (LA-LOO)');

subplot(2,2,2), plot(pred.Eft, Eft_cav,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(cavity) vs. LOO-CV');
xlabel('Latent prediction (CV)'); ylabel('Latent prediction (LA-LOO)');

subplot(2,2,3), plot(pred.Eft, Eft_inla,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('LA-LOO(inla) vs. LOO-CV');
xlabel('Latent prediction (CV)'); ylabel('Latent prediction (LA-LOO)');

subplot(2,2,4), plot(pred_ep.Eft, Eft_ep,'.'); ax = axis;
line([min(ax(1),ax(3)) max(ax(2),ax(4))], [min(ax(1),ax(3)) max(ax(2),ax(4))]);
title('EP-LOO vs. LOO-CV');
xlabel('Latent prediction (CV)'); ylabel('Latent prediction (LA-LOO)');
