%DEMO_CLASSIFIC  Classification problem demonstration for 2 classes
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
%    Here we demonstarte use of Laplace, EP and MCMC methods to
%    find the posterior of the latent values and parameters. With
%    these we can make predictions on the class probability of
%    future observations. See Rasmussen & Willimas (2006) for the
%    detailed treatment of Laplace and EP for probit and logit
%    models and Neal (1998) for MCMC approach for probit and logit
%    models.
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

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

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

% Make predictions
[Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = ...
    gp_pred(gp, x, y, xt, 'yt', ones(size(xt,1),1) );

% Plot some nice figures that show results

% Visualise predictive probability p(ystar = 1) with grayscale
figure, hold on;
n_pred=size(xt,1);
h1=pcolor(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpyt_la),20,20));
set(h1, 'edgealpha', 0), set(h1, 'facecolor', 'interp')
colormap(repmat(linspace(1,0,64)', 1, 3).*repmat(ones(1,3), 64,1))
axis([-inf inf -inf inf]), %axis off
plot(x(y==-1,1),x(y==-1,2),'o', 'markersize', 8, 'linewidth', 2);
plot(x(y==1,1),x(y==1,2),'rx', 'markersize', 8, 'linewidth', 2);
set(gcf, 'color', 'w'), title('predictive probability and training cases with Laplace', 'fontsize', 14)

% Visualise predictive probability p(ystar = 1) with contours
figure, hold on
[cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpyt_la),20,20),[0.025 0.25 0.5 0.75 0.975], 'linewidth', 3);
text_handle = clabel(cs,h);
set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7],'linewidth', 2, 'fontsize',14)
c1=[linspace(0,1,64)' 0*ones(64,1) linspace(1,0,64)'];
colormap(c1)
plot(x(y==1,1), x(y==1,2), 'rx', 'markersize', 8, 'linewidth', 2),
plot(x(y==-1,1), x(y==-1,2), 'bo', 'markersize', 8, 'linewidth', 2)
plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
set(gcf, 'color', 'w'), title('predictive probability contours with Laplace', 'fontsize', 14)


% ------- Expectation propagation --------
fprintf(['%s model with EP integration over the latent values and\n' ...
         'MAP estimate for the parameters\n'],gp.lik.type)

% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Make predictions
[Eft_ep, Varft_ep, lpyt_ep, Eyt_ep, Varyt_ep] = ...
    gp_pred(gp, x, y, xt, 'yt', ones(size(xt,1),1) );

% Plot some nice figures that show results

% Visualise predictive probability p(ystar = 1) with grayscale
figure, hold on;
n_pred=size(xt,1);
h1=pcolor(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpyt_ep),20,20));
set(h1, 'edgealpha', 0), set(h1, 'facecolor', 'interp')
colormap(repmat(linspace(1,0,64)', 1, 3).*repmat(ones(1,3), 64,1))
axis([-inf inf -inf inf]), %axis off
plot(x(y==-1,1),x(y==-1,2),'o', 'markersize', 8, 'linewidth', 2);
plot(x(y==1,1),x(y==1,2),'rx', 'markersize', 8, 'linewidth', 2);
set(gcf, 'color', 'w'), title('predictive probability and training cases with EP', 'fontsize', 14)

% Visualise predictive probability  p(ystar = 1) with contours
figure, hold on
[cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpyt_ep),20,20),[0.025 0.25 0.5 0.75 0.975], 'linewidth', 3);
text_handle = clabel(cs,h);
set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7],'linewidth', 2, 'fontsize',14)
c1=[linspace(0,1,64)' 0*ones(64,1) linspace(1,0,64)'];
colormap(c1)
plot(x(y==1,1), x(y==1,2), 'rx', 'markersize', 8, 'linewidth', 2),
plot(x(y==-1,1), x(y==-1,2), 'bo', 'markersize', 8, 'linewidth', 2)
plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
set(gcf, 'color', 'w'), title('predictive probability contours with EP', 'fontsize', 14)


% ------- MCMC ---------------
fprintf(['%s model with MCMC integration over the latent values and\n' ...
         'the parameters\n'],gp.lik.type)

% Set the approximate inference method
% Note that MCMC for latent values requires often more jitter
gp = gp_set(gp, 'latent_method', 'MCMC', 'jitterSigma2', 1e-6);

% Sample using default method, that is, surrogate and elliptical slice samplers
% these samplers are quite robust with default options
[gp_rec,g,opt]=gp_mc(gp, x, y, 'nsamples', 220, 'display', 20);
% Remove burn-in and thin
gp_rec=thin(gp_rec,21,2);

% Make predictions
[Ef_mc, Varf_mc, lpy_mc, Ey_mc, Vary_mc] = ...
    gp_pred(gp_rec, x, y, xt, 'yt', ones(size(xt,1),1) );

% Plot some nice figures that show results

% Visualise predictive probability p(ystar = 1) with grayscale
figure, hold on;
n_pred=size(xt,1);
h1=pcolor(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpy_mc),20,20));
set(h1, 'edgealpha', 0), set(h1, 'facecolor', 'interp')
colormap(repmat(linspace(1,0,64)', 1, 3).*repmat(ones(1,3), 64,1))
axis([-inf inf -inf inf]), %axis off
plot(x(y==-1,1),x(y==-1,2),'o', 'markersize', 8, 'linewidth', 2);
plot(x(y==1,1),x(y==1,2),'rx', 'markersize', 8, 'linewidth', 2);
set(gcf, 'color', 'w'), title('predictive probability and training cases with MCMC', 'fontsize', 14)

% Visualise predictive probability  p(ystar = 1) with contours
figure, hold on
[cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(exp(lpy_mc),20,20),[0.025 0.25 0.5 0.75 0.975], 'linewidth', 3);
text_handle = clabel(cs,h);
set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7],'linewidth', 2, 'fontsize',14)
c1=[linspace(0,1,64)' 0*ones(64,1) linspace(1,0,64)'];
colormap(c1)
plot(x(y==1,1), x(y==1,2), 'rx', 'markersize', 8, 'linewidth', 2),
plot(x(y==-1,1), x(y==-1,2), 'bo', 'markersize', 8, 'linewidth', 2)
plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
set(gcf, 'color', 'w'), title('predictive probability contours with MCMC', 'fontsize', 14)

% ------- Comparison ---------------
disp('Compare MCMC, Laplace and EP results for two latent variables')

% compare MCMC, Laplace and EP results for two latent variables
% here MCMC result includes uncertainty related to parameters,
% while Laplace and EP results use MAP value for parameters
% see GP_IA for integrating over parameters when using Laplace
% or EP for latent values
apu1 = 123; apu2 = 340;
[Efs_mc, Varfs_mc, lpys_mc] = gpmc_preds(gp_rec, x, y, xt, 'yt', ones(size(xt,1),1) );
sf1 = randn(size(Efs_mc(apu1,:))).*sqrt(Varfs_mc(apu1,:))+Efs_mc(apu1,:);
sf2 = randn(size(Efs_mc(apu2,:))).*sqrt(Varfs_mc(apu2,:))+Efs_mc(apu2,:);

figure

subplot(1,2,1)
[N,X] = hist(sf1,20);
hist(sf1,20);
h1=get(gca,'Children');
hold on
x_in = min(sf1)-2:0.1:max(sf1)+4;
ff = norm_pdf(x_in, Eft_la(apu1), sqrt(Varft_la(apu1)));
h2=plot(x_in, max(N)/max(ff)*ff, 'g', 'lineWidth', 2);
ff = norm_pdf(x_in, Eft_ep(apu1), sqrt(Varft_ep(apu1)));
h3=plot(x_in, max(N)/max(ff)*ff, 'r', 'lineWidth', 2);
set(gca, 'Ytick', [])
legend([h1 h2 h3],'MCMC','Laplace','EP')
title(sprintf('p(f|D) at input location (%.1f, %.1f)', xt(apu1,1), xt(apu1,2)));
xlim([-15 5])

subplot(1,2,2)
[N,X] = hist(sf2,20);
hist(sf2,20)
h1=get(gca,'Children');
hold on
x_in = min(sf2)-2:0.1:max(sf2)+2;
ff = norm_pdf(x_in, Eft_la(apu2), sqrt(Varft_la(apu2)));
h2=plot(x_in, max(N)/max(ff)*ff, 'g', 'lineWidth', 2);
ff = norm_pdf(x_in, Eft_ep(apu2), sqrt(Varft_ep(apu2)));
h3=plot(x_in, max(N)/max(ff)*ff, 'r', 'lineWidth', 2);
set(gca, 'Ytick', [])
legend([h1 h2 h3],'MCMC','Laplace','EP')
title(sprintf('p(f|D) at input location (%.1f, %.1f)', xt(apu2,1), xt(apu2,2)));
xlim([-2 10])







% $$$ % ======================
% $$$ % Print the figures for the manual
% $$$ % ======================
% $$$ figure, hold on
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_la,20,20),[0.75 0.975], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_la,20,20),[0.5], 'linewidth', 2.5, 'color', 'k', 'lineStyle', '-');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_la,20,20),[0.025 0.25], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ plot(x(y==1,1), x(y==1,2), 'kx', 'markersize', 4, 'linewidth', 1),
% $$$ plot(x(y==-1,1), x(y==-1,2), 'ko', 'markersize', 2, 'linewidth', 1)
% $$$ %plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
% $$$ set(gcf, 'color', 'w')
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 5 4])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_classific1_figLA.eps
% $$$ 
% $$$ figure, hold on
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_ep,20,20),[0.75 0.975], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_ep,20,20),[0.5], 'linewidth', 2.5, 'color', 'k', 'lineStyle', '-');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Py_ep,20,20),[0.025 0.25], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ plot(x(y==1,1), x(y==1,2), 'kx', 'markersize', 4, 'linewidth', 1),
% $$$ plot(x(y==-1,1), x(y==-1,2), 'ko', 'markersize', 2, 'linewidth', 1)
% $$$ %plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
% $$$ set(gcf, 'color', 'w')
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 5 4])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_classific1_figEP.eps
% $$$ 
% $$$ figure, hold on
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Pyt_mc,20,20),[0.75 0.975], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Pyt_mc,20,20),[0.5], 'linewidth', 2.5, 'color', 'k', 'lineStyle', '-');
% $$$ [cs,h]=contour(reshape(xt(:,1),20,20),reshape(xt(:,2),20,20),reshape(Pyt_mc,20,20),[0.025 0.25], 'linewidth', 1, 'color', 'k', 'lineStyle', '--');
% $$$ plot(x(y==1,1), x(y==1,2), 'kx', 'markersize', 4, 'linewidth', 1),
% $$$ plot(x(y==-1,1), x(y==-1,2), 'ko', 'markersize', 2, 'linewidth', 1)
% $$$ %plot(xt(:,1), xt(:,2), 'k.'), axis([-inf inf -inf inf]), %axis off
% $$$ set(gcf, 'color', 'w')
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 5 4])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_classific1_figMCMC.eps
% $$$ 
% $$$ 
% $$$ 
% $$$ % compare MCMC, Laplace and EP results for two latent variables
% $$$ apu1 = 123; apu2 = 340;
% $$$ %apu1 = randpick(1:400);  apu2 = randpick(1:400);
% $$$ sf = Ef_mc(apu1,:);
% $$$ sf2 = Ef_mc(apu2,:);
% $$$ 
% $$$ figure
% $$$ subplot(1,2,1)
% $$$ [N,X] = hist(sf);
% $$$ hist(sf)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ x_in = min(sf)-2:0.1:max(sf)+4;
% $$$ ff = norm_pdf(x_in, Ef_la(apu1), sqrt(Varf_la(apu1)));
% $$$ plot(x_in, max(N)/max(ff)*ff, 'k--', 'lineWidth', 2)
% $$$ ff = norm_pdf(x_in, Ef_ep(apu1), sqrt(Varf_ep(apu1)));
% $$$ plot(x_in, max(N)/max(ff)*ff, 'k', 'lineWidth', 2)
% $$$ %ylim([0 105])
% $$$ set(gca, 'Ytick', [])
% $$$ xlim([-15 5])
% $$$ 
% $$$ subplot(1,2,2)
% $$$ [N,X] = hist(sf2);
% $$$ hist(sf2)
% $$$ h = findobj(gca,'Type','patch');
% $$$ set(h,'FaceColor','w','EdgeColor','k')
% $$$ hold on
% $$$ x_in = min(sf2)-2:0.1:max(sf2)+2;
% $$$ ff = norm_pdf(x_in, Ef_la(apu2), sqrt(Varf_la(apu2)));
% $$$ plot(x_in, max(N)/max(ff)*ff, 'k--', 'lineWidth', 2)
% $$$ ff = norm_pdf(x_in, Ef_ep(apu2), sqrt(Varf_ep(apu2)));
% $$$ plot(x_in, max(N)/max(ff)*ff, 'k', 'lineWidth', 2)
% $$$ %ylim([0 105])
% $$$ set(gca, 'Ytick', [])
% $$$ xlim([-2 10])
% $$$ 
% $$$ 
% $$$ 
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'pos',[15 14 7 5])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ print -depsc2 /proj/bayes/jpvanhat/software/doc/GPstuffDoc/pics/demo_classific1_figHist.eps
% $$$ 

