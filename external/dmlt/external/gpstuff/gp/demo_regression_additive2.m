%DEMO_REGRESSION_ADDITIVE2  Regression demonstration with additive Gaussian
%                           process using linear, squared exponential and
%                           neural network covariance fucntions 
%
%  Description
%    Gaussian process solutions in 2D regression problem using
%    constant, linear, squared exponential (sexp) and neural
%    network covariance functions, and with various additive
%    combinations of these four covariance functions. The noisy
%    observations y are assumed to satisfy
%
%         y = f + e,    where e ~ N(0, s^2)
%
%    where f is an unknown underlying function. A zero mean Gaussian
%    process prior is assumed for f
%
%         f ~ N(0, K),
%
%    where K is the covariance matrix whose elements are given by one of
%    the following six covariance function:
%    
%    - constant + linear
%    - costant + sexp for 1. input + linear for 2. input
%    - sexp for 1. input + sexp for 2. input
%    - sexp
%    - neural network for 1. input + neural network for 2. input
%    - neural network
%
%    A prior is assumed for parameters of the covariance functions,
%    and the inference is done with a MAP estimate for parameter
%    values.
%
%    For more detailed discussion of  covariance functions, see e.g.
%
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%
%  See also
%    DEMO_REGRESSION1

% Copyright (c) 2010 Jaakko Riihimäki, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


% REGRESSION TOY DATA
% Create an example regression data
x=rand(225,2)*2-1;
y=3*norm_cdf(2*x(:,2))+2*norm_cdf(4*x(:,1));
% add some noise
y=y+randn(size(y))*0.25;
y=y-mean(y);
[n, nin] = size(x);
% create equally spaced points to visualise the predictions:
[xt1,xt2]=meshgrid(-2:0.1:2,-2:0.1:2);
xt=[xt1(:) xt2(:)];
nxt=size(xt1,1);

% Assume a Student-t distribution for the GP parameters
pt = prior_t('nu', 4, 's2', 10);	% a prior structure

% Create a Gaussian noise model
lik = lik_gaussian('sigma2', 0.2^2);

% Set a small amount of jitter 
jitter=1e-4;

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);

% CONSTANT + LINEAR COVARIANCE FUNCTION
disp('Constant + linear covariance function')
% constant covariance function
gpcf_c = gpcf_constant('constSigma2', 1, 'constSigma2_prior', pt);
% linear covariance function
gpcf_l = gpcf_linear('coeffSigma2_prior', pt);
gp = gp_set('lik', lik, 'cf', {gpcf_c gpcf_l}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (constant + linear)');

% CONSTANT + SQUARED EXPONENTIAL COVARIANCE FUNCTION (W.R.T. THE
% FIRST INPUT DIMENSION) + LINEAR (W.R.T. THE SECOND INPUT
% DIMENSION)
fprintf(['Constant + squared exponential covariance function\n' ...
         '(w.r.t. the first input dimension) + linear (w.r.t.\n' ...
         'the second input dimension)\n'])

%Covariance function for the first input variable
gpcf_s1 = gpcf_sexp('selectedVariables', 1, 'lengthScale',0.5, ...
                    'lengthScale_prior', pt, 'magnSigma2', 0.15, ...
                    'magnSigma2_prior', pt);
% gpcf_s1 can be construted also as
%metric1 = metric_euclidean('components', {[1]}, 'lengthScale',[0.5], ...
%                           'lengthScale_prior', pt);
%gpcf_s1 = gpcf_sexp('magnSigma2', 0.15, 'magnSigma2_prior', pt, ...
%                    'metric', metric1);
%Covariance function for the second input variable
gpcf_l2 = gpcf_linear('selectedVariables', 2, 'coeffSigma2_prior', pt);
gp = gp_set('lik', lik, 'cf', {gpcf_c gpcf_s1 gpcf_l2}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (sexp for 1. input + linear for 2. input )');

% ADDITIVE SQUARED EXPONENTIAL COVARIANCE FUNCTION
disp('Additive squared exponential covariance function')

% Covariance function for the second input variable
gpcf_s2 = gpcf_sexp('selectedVariables', 2,'lengthScale',0.5, ...
                    'lengthScale_prior', pt, 'magnSigma2', 0.15, ...
                    'magnSigma2_prior', pt);
% gpcf_s2 can be construted also as
%metric2 = metric_euclidean('components', {[2]},'lengthScale',[0.5], 'lengthScale_prior', pt);
%gpcf_s2 = gpcf_sexp('magnSigma2', 0.15, 'magnSigma2_prior', pt, 'metric', metric2);
gp = gp_set('lik', lik, 'cf', {gpcf_s1,gpcf_s2}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (additive sexp)');

% SQUARED EXPONENTIAL COVARIANCE FUNCTION
disp('Squared exponential covariance function')

gpcf_s = gpcf_sexp('lengthScale', ones(1,nin), 'magnSigma2', 0.2^2, ...
                   'lengthScale_prior', pt, 'magnSigma2_prior', pt);
gp = gp_set('lik', lik, 'cf', {gpcf_s}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (sexp)');

% ADDITIVE NEURAL NETWORK COVARIANCE FUNCTION
disp('Additive neural network covariance function');

gpcf_nn1 = gpcf_neuralnetwork('weightSigma2', 1, 'biasSigma2', 1, 'selectedVariables', [1], ...
                              'weightSigma2_prior', pt, 'biasSigma2_prior', pt);
gpcf_nn2 = gpcf_neuralnetwork('weightSigma2', 1, 'biasSigma2', 1, 'selectedVariables', [2], ...
                              'weightSigma2_prior', pt, 'biasSigma2_prior', pt);
gp = gp_set('lik', lik, 'cf', {gpcf_nn1,gpcf_nn2}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (additive neural network)');

% NEURAL NETWORK COVARIANCE FUNCTION
disp('Neural network covariance function')

gpcf_nn = gpcf_neuralnetwork('weightSigma2', ones(1,nin), 'biasSigma2', 1, ...
                             'weightSigma2_prior', pt, 'biasSigma2_prior', pt);
gp = gp_set('lik', lik, 'cf', {gpcf_nn}, 'jitterSigma2', jitter);

% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Compute predictions in a grid using the MAP estimate
Eft_map = gp_pred(gp, x, y, xt);

% Plot the prediction and data
figure
set(gcf, 'color', 'w')
mesh(xt1, xt2, reshape(Eft_map,nxt,nxt));
hold on
plot3(x(:,1), x(:,2), y, '*')
xlabel('x_1'), ylabel('x_2')
title('The predicted underlying function (neural network)');
