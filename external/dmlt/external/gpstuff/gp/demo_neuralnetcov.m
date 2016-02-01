%DEMO_NEURALNETCOV  Demonstration of Gaussian process with a neural
%                   network covariance function
%                    
%  Description
%    Infinite neural network solutions in 2D and 1D regression
%    problems with a comparison to Gaussian process solution given
%    by squared exponential covariance function. The noisy
%    observations y are assumed to satisfy
%
%         y = f + e,    where e ~ N(0, s^2)
%
%    where f is an unknown underlying function. A zero mean
%    Gaussian process prior is assumed for f
%
%         f ~ N(0, K),
%
%    where K is the covariance matrix whose elements are given by
%    neural network (or squared exponential) covariance function. A
%    prior is assumed for parameters of the covariance functions,
%    and the inference is done with a MAP estimate for parameter
%    values.
%
%    For more detailed discussion of infinite neural networks, see
%    e.g.
%
%      Neal, R. M. (1996). Bayesian Learning for Neural Networks. 
%      Springer-Verlag.
%
%      Williams, C. K. I. (1996). Computing with infinite networks. 
%      In Advances in Neural Information Processing Systems 9. MIT
%      Press, Cambridge, MA.
%
%
%  See also DEMO_REGRESSION1

% Copyright (c) 2010 Jaakko Riihim�ki, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% 2D REGRESSION DATA

% create 2D example data
x=rand(300,2)*2-1;
y=zeros(size(x,1),1); y(x(:,1)>0&x(:,2)>0)=1;
y=y+0.1*randn(size(y));

[n, nin] = size(x);

% --- Construct the model ---

% squared exponential covariance function
gpcf1 = gpcf_sexp('lengthScale', ones(1,nin), 'magnSigma2', 1);
% neural network covariance function
gpcf2 = gpcf_neuralnetwork('weightSigma2', ones(1,nin), 'biasSigma2', 1);
% Gaussian noise structures
lik = lik_gaussian('sigma2', 0.2^2);

% a prior structure for GP parameters
pt = prior_t('s2', 4);
gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pt, 'magnSigma2_prior', pt);
gpcf2 = gpcf_neuralnetwork(gpcf2, 'weightSigma2_prior', pt, 'biasSigma2_prior', pt);

gp = gp_set('lik', lik, 'cf', gpcf1);
gp2 = gp_set('lik', lik, 'cf', gpcf2);

% --- MAP estimate using scaled conjugate gradient algorithm ---
%     (see scg for more details)

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);
gp2=gp_optim(gp2,x,y,'opt',opt);

% create points where predictions are made
[xt1,xt2]=meshgrid(-1.5:0.05:1.5,-1.5:0.05:1.5);
xt=[xt1(:) xt2(:)];
% compute the predictions
[Eft_map, Varft_map] = gp_pred(gp, x, y, xt);
[Eft_map2, Varft_map2] = gp_pred(gp2, x, y, xt);

% Plot the predictions and data
figure, set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,size(xt1,1),size(xt1,2)));
hold on
plot3(x(:,1), x(:,2), y, '*')
axis on;
title('GP (squared exponential) predictions and the data points');

figure, set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map2,size(xt1,1),size(xt1,2)));
hold on
plot3(x(:,1), x(:,2), y, '*')
axis on;
title('GP (neural network) predictions and the data points');

% 1D REGRESSION DATA

% create a 1D toy data
x=rand(100,1)*4-2;
y=norm_pdf(4*x)+0.05*randn(size(x));
[n, nin] = size(x);

gpcf1 = gpcf_sexp('lengthScale', ones(1,nin), 'magnSigma2', 1);
gpcf2 = gpcf_neuralnetwork('weightSigma2', ones(1,nin), 'biasSigma2', 1);

gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pt, 'magnSigma2_prior', pt);
gpcf2 = gpcf_neuralnetwork(gpcf2, 'weightSigma2_prior', pt, 'biasSigma2_prior', pt);
lik = lik_gaussian();

gp = gp_set('lik', lik, 'cf', gpcf1);
gp2 = gp_set('lik', lik, 'cf', gpcf2);

% --- MAP estimate using scaled conjugate gradient algorithm ---
%     (see fminscg for more details)

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);
gp2=gp_optim(gp2,x,y,'opt',opt);

% create points where predictions are made
xgrid=linspace(min(x)-1.5,max(x)+1.5,200)';
[Eft_map, Varft_map, lpyt_map, Eyt_map, Varyt_map] = gp_pred(gp, x, y, xgrid, 'yt', ones(200,1));
[Eft_map2, Varft_map2, lpyt_map2, Eyt_map2, Varyt_map2] = gp_pred(gp2, x, y, xgrid, 'yt', ones(200,1));

% Plot the predictions and data
color1=ones(1,3)*0.8; color2=ones(1,3)*0.5;
figure, set(gcf, 'color', 'w'), hold on
h1=fill([xgrid' fliplr(xgrid')], [(Eyt_map+1.96*sqrt(Varyt_map))' fliplr((Eyt_map-1.96*sqrt(Varyt_map))')], color1, 'edgecolor', color1);
% GP mean
h2=plot(xgrid, Eyt_map, 'color', color2, 'linewidth', 3);
% observations
h3=plot(x, y, 'xk', 'markersize', 10, 'linewidth', 2);
% true function
h4=plot(xgrid, norm_pdf(4*xgrid), 'color', 'r', 'linewidth', 2);
legend([h1 h2 h3 h4], 'GP 95% CI', 'GP mean', 'observations', 'true latent function')
title('GP (squared exponential) predictions and the data points');

figure, set(gcf, 'color', 'w'), hold on
h1=fill([xgrid' fliplr(xgrid')], [(Eyt_map2+1.96*sqrt(Varyt_map2))' fliplr((Eyt_map2-1.96*sqrt(Varyt_map2))')], color1, 'edgecolor', color1);
% GP mean
h2=plot(xgrid, Eyt_map2, 'color', color2, 'linewidth', 3);
% observations
h3=plot(x, y, 'xk', 'markersize', 10, 'linewidth', 2);
% true function
h4=plot(xgrid, norm_pdf(4*xgrid), 'color', 'r', 'linewidth', 2);
legend([h1 h2 h3 h4], 'GP 95% CI', 'GP mean', 'observations', 'true latent function')
title('GP (neural network) predictions and the data points');
