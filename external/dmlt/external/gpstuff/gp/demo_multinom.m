%DEMO_MULTINOM   Demonstration of Gaussian process model with multinomial
%                likelihood with 3 classes
%
%  Description
%    Demonstration of estimating the unknown population proportions
%    in multinomial model with three classes from a spatially distributed
%    trials. Data consists of observations Y describing the counts in each
%    class at spatial locations X. The multinomial model is 
%
%      Y_i ~ Multinomial(Y_i | sum(Y_i), p_i),
%
%    where the vector p_i collects the proportion of success for each
%    class. The total number of trials at location X_i is sum(Y_i). 
%
%    The latent values for N training points and C classes are
%    f=(f1_1,f2_1,...,fN_1,f1_2,f2_2,...,fN_2,...,f1_C,f2_C,...,fN_C)^T,
%    and are given a zero mean Gaussian process prior
%      
%      f ~ N(0, K),
%
%    where K is a block diagonal covariance matrix with blocks
%    K_1,...,K_C whose elements are given by K_ij = k(x_i, x_j |
%    th). The function k(x_i, x_j | th) is covariance function and
%    th its parameters.
%
%    In this demo we approximate the posterior distribution with
%    Laplace approximation and MCMC.

% Reference:
% Teppo Juntunen, Jarno Vanhatalo, Heikki Peltonen and Samu Mäntyniemi
% (2012). Bayesian spatial multispecies modelling to assess pelagic fish
% stocks from acoustic- and trawl-survey data. ICES Journal of Marine
% Science, 69: 95-104.

% Copyright (c) 2010-2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Simulate data 
% ---------------------------------------------------------------------
% Each class has different prior for its latent function
[x1,x2] = meshgrid(1:0.5:10,1:0.5:10);
xt = [x1(:), x2(:)];

gpcf1 = gpcf_sexp('lengthScale', 2, 'magnSigma2', 2);
gpcf2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 3);
gpcf3 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 1);
gp1 = gp_set('lik', lik_gaussian, 'cf', gpcf1, 'jitterSigma2', 1e-4);
gp2 = gp_set('lik', lik_gaussian, 'cf', gpcf2, 'jitterSigma2', 1e-4);
gp3 = gp_set('lik', lik_gaussian, 'cf', gpcf3, 'jitterSigma2', 1e-4);
K1 = gp_trcov(gp1,xt);
K2 = gp_trcov(gp2,xt);
K3 = gp_trcov(gp3,xt);

L1 = chol(K1,'lower');L2 = chol(K2,'lower');L3 = chol(K3,'lower');
f = [L1*randn(length(xt), 1) L2*randn(length(xt), 1) L3*randn(length(xt), 1)];
expf = exp(f);

p = expf./ repmat(sum(expf,2),1,size(expf,2));
n = floor(1+ 1000*rand(size(p,1),1));
yt = mnrnd(n,p);

I = floor((size(xt,1)-1)*rand(80,1)+1);
I = unique(I,'rows');
y = yt(I,:);       % Training y
x = xt(I,:);       % Training x

% Plot the real surface
figure(1), set(gcf, 'color', 'w'), hold on, mc1=mapcolor(f(:,1));mc2=mapcolor(f(:,2));mc3=mapcolor(f(:,3));
subplot(2,3,1);pcolor(x1, x2, reshape(f(:,1),size(x1))),shading flat,colormap(mc1),colorbar; title('real latent 1')
subplot(2,3,2);pcolor(x1, x2, reshape(f(:,2),size(x1))),shading flat,colormap(mc2),colorbar; title('real latent 2')
subplot(2,3,3);pcolor(x1, x2, reshape(f(:,3),size(x1))),shading flat,colormap(mc3),colorbar; title('real latent 3')

% %Plot with contours
% figure(1), set(gcf, 'color', 'w'), hold on
% subplot(2,3,1);contour(x1, x2, reshape(f(:,1),size(x1)),'r', 'linewidth', 2); title('real latent 1')
% subplot(2,3,2);contour(x1, x2, reshape(f(:,2),size(x1)),'b', 'linewidth', 2); title('real latent 2')
% subplot(2,3,3);contour(x1, x2, reshape(f(:,3),size(x1)),'k', 'linewidth', 2); title('real latent 3')
% ---------------------------------------------------------------------

% Inference with Laplace approximation
% ====================================

% Create the model to infer the above data
lik = lik_multinom;
gpcf1 = gpcf_sexp('lengthScale', 2, 'magnSigma2', 2, 'lengthScale_prior', prior_t('nu',4, 's2', 3));
gpcf2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 3, 'lengthScale_prior', prior_t('nu',4, 's2', 3));
gpcf3 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 1, 'lengthScale_prior', prior_t('nu',4, 's2', 3));
gp = gp_set('lik', lik_multinom, 'cf', {gpcf1 gpcf2 gpcf3}, 'jitterSigma2', 1e-4, 'comp_cf', {1 2 3});

% NOTE! if multiple covariance functions per output is used define
% gp.comp_cf, for example, as follows:
% gp.comp_cf = {[1 2] [3 4] [5 6]};
% To see how the model looks like now type gp in the command prompt

% Optimize with the scaled conjugate gradient method
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Display','iter','MaxIter',100);
gp=gp_optim(gp,x,y,'opt',opt);

% make the prediction for test points
[Eft] = gp_pred(gp, x, y, xt, 'yt', ones(size(yt)));
Eft = reshape(Eft, size(xt,1), size(yt,2));

% Plot the result
figure(1)
subplot(2,3,4);pcolor(x1, x2, reshape(Eft(:,1),size(x1))),shading flat,colormap(mc1),colorbar
hold on; plot(x(:,1),x(:,2),'k.'); title('Model pred. (Laplace)')
subplot(2,3,5);pcolor(x1, x2, reshape(Eft(:,2),size(x1))),shading flat,colormap(mc2),colorbar
hold on; plot(x(:,1),x(:,2),'k.'); title('Model pred. (Laplace)')
subplot(2,3,6);pcolor(x1, x2, reshape(Eft(:,3),size(x1))),shading flat,colormap(mc3),colorbar
hold on; plot(x(:,1),x(:,2),'k.'); title('Model pred. (Laplace)')

% %Plot with contours
% figure(2)
% subplot(2,3,4);contour(x1, x2, reshape(Eft(:,1),size(x1)),'r', 'linewidth', 2)
% hold on; plot(x(:,1),x(:,2),'.'); title('Model pred. (Laplace)')
% subplot(2,3,5);contour(x1, x2, reshape(Eft(:,2),size(x1)),'b', 'linewidth', 2)
% hold on; plot(x(:,1),x(:,2),'.'); title('Model pred. (Laplace)')
% subplot(2,3,6);contour(x1, x2, reshape(Eft(:,3),size(x1)),'k', 'linewidth', 2)
% hold on; plot(x(:,1),x(:,2),'.'); title('Model pred. (Laplace)')

% Plot the relative abundances
pyt2=exp(Eft)./(sum(exp(Eft),2)*ones(1,3));

figure, set(gcf, 'color', 'w'), hold on, mc1=mapcolor(p(:,1));mc2=mapcolor(p(:,2));mc3=mapcolor(p(:,3));
subplot(2,3,1);pcolor(x1, x2, reshape(p(:,1),size(x1))),shading flat,colormap(mc1),colorbar
title('Real relative abundance 1')
subplot(2,3,2);pcolor(x1, x2, reshape(p(:,2),size(x1))),shading flat,colormap(mc2),colorbar
title('Real relative abundance 2')
subplot(2,3,3);pcolor(x1, x2, reshape(p(:,3),size(x1))),shading flat,colormap(mc3),colorbar
title('Real relative abundance 3')

subplot(2,3,4);pcolor(x1, x2, reshape(pyt2(:,1),size(x1))),shading flat,colormap(mc1),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model pred. (Laplace)')
subplot(2,3,5);pcolor(x1, x2, reshape(pyt2(:,2),size(x1))),shading flat,colormap(mc2),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model pred. (Laplace)')
subplot(2,3,6);pcolor(x1, x2, reshape(pyt2(:,3),size(x1))),shading flat,colormap(mc3),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model pred. (Laplace)')

%% Inference with MCMC
% ===================================================================

% Set the approximate inference method to MCMC
% Note that MCMC for latent values requires often more jitter
lat = gp_pred(gp, x, y, x);
gp = gp_set(gp, 'latent_method', 'MCMC', 'jitterSigma2', 1e-4);
gp = gp_set(gp, 'latent_opt', struct('method',@scaled_mh));
gp.latentValues = lat(:);      % initialize the latent values

% Set the options for MCMC...
hmc_opt.steps=10;
hmc_opt.stepadj=0.001;
hmc_opt.nsamples=1;
latent_opt.display=0;
latent_opt.repeat = 20;
latent_opt.sample_latent_scale = 0.05;
hmc2('state', sum(100*clock))

% Sample an initial position
[r,g,opt]=gp_mc(gp, x, y, 'hmc_opt', hmc_opt, 'latent_opt', latent_opt, 'nsamples', 1, 'repeat', 15);

% re-set some of the sampling options to be more efficient
hmc_opt.repeat=1;
hmc_opt.steps=4;
hmc_opt.stepadj=0.03;
latent_opt.repeat = 5;
hmc2('state', sum(100*clock));

% Sample 
[rgp,g,opt]=gp_mc(g, x, y, 'nsamples', 400, 'hmc_opt', hmc_opt, 'latent_opt', latent_opt, 'record', r);
% Remove burn-in
rgp=thin(rgp,102,2);

% Make predictions
Efs_mc = gpmc_preds(rgp, x, y, xt);
% [Efs_mc, Varfs_mc, lpgs_mc] = gpmc_mo_preds(rgp, x, y, xt, 'yt', ones(size(xt,1),3));

Ef_mc = reshape(mean(Efs_mc,2),361,3);
% pg_mc = reshape(mean(exp(lpgs_mc),2),361,3);


pyt2_mc = exp(Ef_mc)./(sum(exp(Ef_mc),2)*ones(1,3));
% pyt2_mc = reshape(mean(pgs_mc,2),900,3);

% Plot the relative abundances
figure, set(gcf, 'color', 'w'), hold on, mc1=mapcolor(p(:,1));mc2=mapcolor(p(:,2));mc3=mapcolor(p(:,3));
subplot(2,3,1);pcolor(x1, x2, reshape(p(:,1),size(x1))),shading flat,colormap(mc1),colorbar
title('Real relative abundance 1')
subplot(2,3,2);pcolor(x1, x2, reshape(p(:,2),size(x1))),shading flat,colormap(mc2),colorbar
title('Real relative abundance 2')
subplot(2,3,3);pcolor(x1, x2, reshape(p(:,3),size(x1))),shading flat,colormap(mc3),colorbar
title('Real relative abundance 3')

subplot(2,3,4);pcolor(x1, x2, reshape(pyt2_mc(:,1),size(x1))),shading flat,colormap(mc1),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model prediction (MCMC)')
subplot(2,3,5);pcolor(x1, x2, reshape(pyt2_mc(:,2),size(x1))),shading flat,colormap(mc2),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model prediction (MCMC)')
subplot(2,3,6);pcolor(x1, x2, reshape(pyt2_mc(:,3),size(x1))),shading flat,colormap(mc3),colorbar
hold on; plot(x(:,1),x(:,2),'k.'), title('Model prediction (MCMC)')
