%DEMO_QUANTILEGP Demonstration of Quantile GP regression
%
%  Description 
%    This demo demonstrates Quantile-GP regression with toy data. QGP can
%    be used for estimating quantiles of interest of the response variable
%    with respect to input variables.
%
%  See also  LIK_QGP, DEMO_*

% Copyright (c) 2013 Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

prevstream=setrandstream(0);
% Toy data
n=200;
x=linspace(-3,2,n)';
y=((x.^3 + 3.*x.^2-6.*x-8)/4);
x=(x-mean(x))/std(x);
y=(y-mean(y))/std(y);
xt=x;
yt=y;

% add heteroscedastic noise
noise=0.2.*randn(n,1);
noise=noise + 0.2.*exp(-1.5.*x).*randn(n,1);
y=y+noise;

% Create covariance function
gpcf=gpcf_sexp('lengthScale', 1, 'magnSigma2', 1);

% Create quantile-GP likelihood 
lik=lik_qgp('sigma2', 0.5, 'sigma2_prior', prior_t('s2',10,'nu',1));

% We want to infer the 0.25, 0.5 and 0.75 quantiles
lik=lik_qgp(lik, 'quantile', 0.05);
lik2=lik_qgp(lik, 'quantile', 0.5);
lik3=lik_qgp(lik, 'quantile', 0.95);

lik4=lik_qgp(lik, 'quantile', 0.25);
lik5=lik_qgp(lik, 'quantile', 0.75);

% Create gp structure. Quantile-GP only works with EP or MCMC. Lets use
% MCMC for 5%, 50% and 95% quantiles.
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-6, 'latent_method', 'MCMC');
gp2 = gp_set('lik', lik2, 'cf', gpcf, 'jitterSigma2', 1e-6, 'latent_method', 'MCMC');
gp3 = gp_set('lik', lik3, 'cf', gpcf, 'jitterSigma2', 1e-6, 'latent_method', 'MCMC');

% EP for 25% and 75% quantiles (with low number of data points, EP is 
% unstable for very low or high quantiles)
gp4 = gp_set('lik', lik4, 'cf', gpcf, 'jitterSigma2', 1e-6, 'latent_method', 'EP');
gp5 = gp_set('lik', lik5, 'cf', gpcf, 'jitterSigma2', 1e-6, 'latent_method', 'EP');

% Set options for optimization
opt=optimset('TolX',1e-4, 'TolFun', 1e-4, 'Display', 'iter', 'derivativecheck', 'off');

% Sample values
ssls_opt.latent_opt.repeat=100;
rgp=gp_mc(gp, x, y, 'nsamples',200,'ssls_opt', ssls_opt, 'display', 10);
rgp=thin(rgp,10);
rgp2=gp_mc(gp2, x, y, 'nsamples', 200, 'ssls_opt', ssls_opt, 'display', 10);
rgp2=thin(rgp2,10);
rgp3=gp_mc(gp3, x, y, 'nsamples', 200, 'ssls_opt', ssls_opt, 'display', 10);
rgp3=thin(rgp3,10);

% Optimize hyperparameters
gp4=gp_optim(gp4, x, y, 'opt', opt, 'optimf', @fminlbfgs);
gp5=gp_optim(gp5, x, y, 'opt', opt, 'optimf', @fminlbfgs);

% Do predictions for training points
[Ef1,Varf1]=gp_pred(rgp, x, y, xt);
[Ef2,Varf2]=gp_pred(rgp2, x, y, xt);
[Ef3,Varf3]=gp_pred(rgp3, x, y, xt);
[Ef4,Varf4]=gp_pred(gp4, x, y, xt);
[Ef5,Varf5]=gp_pred(gp5, x, y, xt);

% Visualize
figure(1);
plot(xt, Ef1, '-r',xt, Ef4,'-m', xt, Ef2, '-k', xt, Ef3, '-r', xt, y, '.',xt, Ef5,'-m');
legend('5% / 95%', '25% / 75%', '50%');
setrandstream(prevstream);