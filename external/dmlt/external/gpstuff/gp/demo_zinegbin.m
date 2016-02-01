%DEMO_ZINEGBIN  Demonstration of zero-inflated Negative-binomial model
%               using Gaussian process prior
%
%  Description
%    Zero-inflated Negative-binomial model provides a way of modelling
%    count data with excess of zeros. The latent values for N training
%    points are f=(f1_1,f2_1,...,fN_1,f1_2,f2_2,...,fN_2)^T, where latents
%    f_1 are associated with classification process and latents f_2 with
%    count process. Both processes are given a zero mean Gaussian process
%    prior
%
%    f ~ N(0, K),
%
%    where K is a block diagonal covariance matrix with blocks K_1, K_2
%    whose elements are given by K_ij = k(x_i, x_j | th). The function
%    k(x_i, x_j | th) is covariance function and th its parameters.  
%
%    In this demo we approximate the posterior distribution with Laplace
%    approximation. 
%
%    See also  DEMO_SPATIAL2, DEMO_CLASSIFIC1

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari
% Copyright (c) 2011 Jaakko Riihimäki

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% load the data
S = which('demo_spatial1');
data = load(strrep(S,'demo_spatial1.m','demodata/spatial1.txt'));

x = data(:,1:2);
ye = data(:,3);
y = data(:,4);

x0=x;
x=bsxfun(@rdivide,bsxfun(@minus,x0,mean(x0)),std(x0));

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_t('s2',10); 
% pm = prior_sqrtunif();

pl = prior_t('s2',10);
pm = prior_sqrtunif();
cf = gpcf_matern32('lengthScale', 5, 'magnSigma2', 0.05, ...
                   'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
lik = lik_zinegbin('disper_prior', prior_fixed());

% NOTE! if multiple covariance functions per latent is used, define
% gp.comp_cf as follows:
% gp = gp_set(..., 'comp_cf' {[1 2] [5 6]};
gp = gp_set('lik', lik, 'cf', {cf cf}, 'jitterSigma2', 1e-6, 'comp_cf', {[1] [2]});

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt, 'optimf', @fminlbfgs);

% make prediction to the data points
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye);

% Define help parameters for plotting
xii=sub2ind([60 35],x0(:,2),x0(:,1));
[X1,X2]=meshgrid(1:35,1:60);

% Plot the figures
figure
G=repmat(NaN,size(X1));
G(xii)=Ef(1:size(x,1));
pcolor(X1,X2,G),shading flat
colorbar
axis equal
axis([0 35 0 60])
title('Posterior mean of latent (classification process)')

% Plot the figures
figure
G=repmat(NaN,size(X1));
G(xii)=(Ef((size(x,1)+1):end));
pcolor(X1,X2,G),shading flat
colorbar
axis equal
axis([0 35 0 60])
title('Posterior mean of latent (count process)')
