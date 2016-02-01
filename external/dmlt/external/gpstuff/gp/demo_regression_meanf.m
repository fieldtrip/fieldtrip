%DEMO_REGRESSION_MEANF  Regression problem demonstration for GP model with a
%                       mean function
%
%  Description
%    The regression problem consist of a data with one input
%    variable and one output variable with Gaussian noise. The
%    problem is modelled with Full GP model with gaussian
%    likelihood and a specified mean function. The mean function m
%    is a weighted sum of some basis functions h, where
%
%                   m=h'*Beta
%
%    and we have set a gaussian prior for the weights Beta
%
%                   Beta ~ N(b,B)
%
%    Inference is done according to Rasmussen and Williams (2006)
%    p. 27-29.
%
%    In this demonstration the data is from a function:
%
%         y = 2 + x + x^2 + 4*cos(x)*sin(x) + epsilon
%
%    and we define the base functions to be:
%
%       h1(x)=x^2, h2(x)=x, h3(x)=2           
%
%    References:
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.

% Copyright (c) 2010 Tuomas Nikoskinen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Create the data
 x=[-2:0.6:2]';
 res=4*cos(x).*sin(x)+0.4*randn(size(sin(x)));
 y= 2 + x + x.^2 + res;
%---------------
 
gpcf = gpcf_sexp('lengthScale', [0.5], 'magnSigma2', .5);
lik = lik_gaussian('sigma2', 0.4^2);

% Initialize base functions for GP's mean function.
gpmf1 = gpmf_constant('prior_mean',.3,'prior_cov',1);
gpmf2 = gpmf_linear('prior_mean',.3,'prior_cov',1);
gpmf3 = gpmf_squared('prior_mean',.3,'prior_cov',1);

% Initialize gp structure
gp = gp_set('lik', lik, 'cf', gpcf, 'meanf', {gpmf1,gpmf2,gpmf3});

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'DerivativeCheck','on');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Predictions
xt=[-3:0.1:3]';
[Eft, Varft] = gp_pred(gp, x, y, xt);

% PLOT THE DATA

figure
m=plot(xt,Eft,'k','lineWidth',2);
hold on
plot(xt,Eft+2*sqrt(Varft),'k--')
hold on
m95=plot(xt,Eft-2*sqrt(Varft),'k--');
hold on
hav=plot(x,y, 'ro','markerSize',6,'MarkerFaceColor','r');
hold on
h=plot(xt,2+xt+xt.^2+4*cos(xt).*sin(xt),'b--','lineWidth',2);
hold on
mmmean=plot(xt,2+xt+xt.^2,'r--','lineWidth',1);
legend([m m95 h mmmean hav],'prediction','95%','f(x)','mean function','observations');
xlabel('input x')
ylabel('output y')
title('GP regression with a mean function')
