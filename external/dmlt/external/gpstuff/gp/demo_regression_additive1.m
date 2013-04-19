%DEMO_REGRESSION_ADDITIVE1  Regression demonstration with additive model
%
%  Description
%    A regression demonstration with one input variable and one output
%    variable with Gaussian noise. The output is assumed to be
%    realization of two additive functions and Gaussian noise.
%
%    The model constructed is following:
%
%    The observations y are assumed to satisfy
%
%         y = f + g + e,    where e ~ N(0, s^2).
%
%    f and g are underlying latent functions, which we are
%    interested in. We place a zero mean Gaussian process prior for
%    them, which implies that at the observed input locations
%    latent values have prior
%
%         f ~ N(0, Kf) and g ~ N(0,Kg)
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters.
%
%    Since both likelihoods and prior are Gaussian, we obtain a
%    Gaussian marginal likelihood
%
%        p(y|th) = N(0, Kf + Kg + I*s^2).
%    
%    By placing a prior for parameters, p(th), we can find the
%    maximum a posterior (MAP) estimate for them by maximizing
%
%       argmax   log p(y|th) + log p(th).
%         th
%
%    After finding MAP estimate or posterior samples of parameters,
%    we can use them to make predictions for the latent functions. 
%    For example, the posterior predictive distribution of f is:
%
%       p(f | y, th) = N(m, S),
%       m = Kf * (Kf + Kg + s^2I)^(-1) * y
%       S = Kf - Kf * (Kf + Kg + s^2I)^(-1) * Kf
%
%    (We could integrate also over the parameters with, for
%    example, grid integration or MCMC. This is not demonstrated
%    here but it is done exactly similarly as in the
%    demo_regression1.)
%   
%    The demo is organised in four parts:
%     1) data analysis with full GP model
%     2) data analysis with FIC approximation
%     3) data analysis with PIC approximation
%     4) data analysis with CS+FIC model
%
%    For more detailed discussion of Gaussian process regression
%    see Rasmussen and Williams (2006) and for a detailed
%    discussion on sparse additive models see Vanhatalo and Vehtari
%    (2008).
%
% References:
%
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%    Vanhatalo, J. and Vehtari, A. (2008). Modelling local and global
%    phenomena with sparse Gaussian processes. Proceedings of the 24th
%    Conference on Uncertainty in Artificial Intelligence.
%
% See also 
%  DEMO_REGRESSION1, DEMO_SPARSEREGRESION
%

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

%========================================================
% PART 1 data analysis with full GP model
%========================================================

% Load the data
S = which('demo_regression1');
L = strrep(S,'demo_regression1.m','demodata/maunaloa_data.txt');
data=load(L);
y = data(:, 2:13);
y=y';
y=y(:);
x = [1:1:length(y)]';
x = x(y>0);             % Remove contaminated observations
y = y(y>0);
avgy = mean(y);
y = y-avgy;
xt = [0:0.5:565]';

[n,nin] = size(x);
% Now 'x' consist of the inputs and 'y' of the output. 
% 'n' and 'nin' are the number of data points and the 
% dimensionality of 'x' (the number of inputs).

% ---------------------------
% --- Construct the model ---
% 
% First create squared exponential and piecewise polynomial 2
% covariance functions and Gaussian noise structures and set
% priors for their parameters (if SuiteSparse is not
% installed, use gpcf_sexp instead of gpcf_ppcs2)
pl1 = prior_t('s2', 100, 'nu', 10);
pl2 = prior_t('s2', 10, 'nu', 10);
pm1 = prior_sqrtt('s2',300);
pm2 = prior_sqrtt('s2',10);
pn = prior_logunif();
gpcf1 = gpcf_sexp('lengthScale', 100, 'magnSigma2', 200, 'lengthScale_prior', pl1, 'magnSigma2_prior', pm1);
if exist('ldlchol')
  gpcf2 = gpcf_ppcs2('nin', nin, 'lengthScale', 5, 'magnSigma2', 5, 'lengthScale_prior', pl2, 'magnSigma2_prior', pm2);
else
  warning('GPstuff:SuiteSparseMissing',...
  ['SuiteSparse is not properly installed. (in BECS try ''use suitesparse'')\n' ...
   'Using gpcf_sexp (non-compact support) instead of gpcf_ppcs2 (compact support)']);
  gpcf2 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 1, 'lengthScale_prior', pl2, 'magnSigma2_prior', pm);
end
lik = lik_gaussian('sigma2', 0.1, 'sigma2_prior', pn);

% Create the GP structure
gp = gp_set('lik', lik, 'cf', {gpcf1, gpcf2}, 'jitterSigma2', 1e-9) 

% -----------------------------
% --- Conduct the inference ---
%
% --- MAP estimate -----------
% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% Make predictions. Below Ef_full is the predictive mean and Varf_full
% the predictive variance.
[Eft_full, Varft_full, lpyt_full, Eyt_full, Varyt_full] = gp_pred(gp, x, y, xt, 'yt', ones(size(xt)));
[Eft_full1, Varft_full1] = gp_pred(gp, x, y, xt, 'predcf', 1);
[Eft_full2, Varft_full2] = gp_pred(gp, x, y, xt, 'predcf', 2);

% Plot the prediction and data
figure
subplot(2,1,1)
hold on
plot(x,y,'.', 'MarkerSize',7)
plot(xt,Eyt_full,'k', 'LineWidth', 2)
plot(xt,Eyt_full-2.*sqrt(Varyt_full),'g--')
plot(xt,Eyt_full+2.*sqrt(Varyt_full),'g--')
axis tight
caption1 = sprintf('Full GP:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f \n s^2_{noise} = %.2f', gp.cf{1}.lengthScale, gp.cf{1}.magnSigma2, gp.cf{2}.lengthScale, gp.cf{2}.magnSigma2, gp.lik.sigma2);
title(caption1)
legend('Data point', 'predicted mean', '2\sigma error',4)

subplot(2,1,2)
[AX, H1, H2] = plotyy(xt, Eft_full2, xt, Eft_full1);
set(H2,'LineStyle','--')
set(H2, 'LineWidth', 2)
%set(H1, 'Color', 'k')
set(H1,'LineStyle','-')
set(H1, 'LineWidth', 0.8)
title('The long and short term trend')

%========================================================
% PART 2 data analysis with FIC approximation
%========================================================

% Here we conduct the same analysis as in part 1, but this time 
% using FIC approximation. Notice that both covariance components 
% utilize the inducing inputs. This leads to problems since the 
% number of inducing inputs is too small to capture the short term 
% variation. In CS+FIC (later model) the compact support function 
% does not utilize the inducing inputs and for this reason it is
% able to capture also the fast variations.

% Place inducing inputs evenly
Xu = [min(x):24:max(x)+10]';

% Create the FIC GP structure
gp_fic = gp_set('type', 'FIC', 'lik', lik, 'cf', {gpcf1,gpcf2}, 'jitterSigma2', 1e-9, 'X_u', Xu)

% -----------------------------
% --- Conduct the inference ---

% --- MAP estimate using modified Newton algorithm ---

% Now you can choose, if you want to optimize only parameters or 
% optimize simultaneously parameters and inducing inputs. Note that 
% the inducing inputs are not transformed through logarithm when packed

%gp_fic = gp_set(gp_fic, 'infer_params', 'covariance+likelihood+inducing');  % optimize parameters and inducing inputs
gp_fic = gp_set(gp_fic, 'infer_params', 'covariance+likelihood'); % optimize only parameters

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp_fic=gp_optim(gp_fic,x,y,'opt',opt);

% Make the prediction
[Eft_fic, Varft_fic, lpyt_fic, Eyt_fic, Varyt_fic] = gp_pred(gp_fic, x, y, xt, 'yt', ones(size(xt)));

% Plot the solution of FIC
figure
%subplot(4,1,1)
hold on
plot(x,y,'.', 'MarkerSize',7)
plot(xt,Eyt_fic,'k', 'LineWidth', 2)
plot(xt,Eyt_fic-2.*sqrt(Varyt_fic),'g--', 'LineWidth', 2)
plot(gp_fic.X_u, -30, 'rx', 'MarkerSize', 5, 'LineWidth', 2)
plot(xt,Eyt_fic+2.*sqrt(Varyt_fic),'g--', 'LineWidth', 2)
axis tight
caption2 = sprintf('FIC:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f \n s^2_{noise} = %.2f', gp_fic.cf{1}.lengthScale, gp_fic.cf{1}.magnSigma2, gp_fic.cf{2}.lengthScale, gp_fic.cf{2}.magnSigma2, gp_fic.lik.sigma2);
title(caption2)
legend('Data point', 'predicted mean', '2\sigma error', 'inducing input','Location','Northwest')


%========================================================
% PART 3 data analysis with PIC approximation
%========================================================

% set the data points into clusters
edges = linspace(-1,max(x)+1,20);
tot=0; 
for i=1:length(edges)-1
    trindex{i} = find(x>edges(i) & x<edges(i+1));
end
% Create the FIC GP structure
gp_pic = gp_set('type', 'PIC', 'lik', lik, 'cf', {gpcf1, gpcf2}, 'jitterSigma2', 1e-6, 'X_u', Xu)
gp_pic = gp_set(gp_pic, 'tr_index', trindex);

% -----------------------------
% --- Conduct the inference ---

% --- MAP estimate using modified Newton algorithm ---

% Now you can choose, if you want to optimize only parameters or 
% optimize simultaneously parameters and inducing inputs. Note that 
% the inducing inputs are not transformed through logarithm when packed

%gp_pic = gp_set(gp_pic, 'infer_params', 'covariance+likelihood+inducing');  % optimize parameters and inducing inputs
gp_pic = gp_set(gp_pic, 'infer_params', 'covariance+likelihood');           % optimize only parameters

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp_pic=gp_optim(gp_pic,x,y,'opt',opt);
gp_pic=gp_optim(gp_pic,x,y,'opt',opt,'optimf',@fminscg);

% Make the prediction
[Eft_pic, Varft_pic, lpyt_pic, Eyt_pic, Varyt_pic] = gp_pred(gp_pic, x, y, x, 'tstind', trindex, 'yt', y);


% Plot the solution of PIC
figure
%subplot(4,1,1)
hold on
plot(x,y,'.', 'MarkerSize',7)
plot(x,Eft_pic,'k', 'LineWidth', 2)
plot(x,Eft_pic-2.*sqrt(Varyt_pic),'g--', 'LineWidth', 2)
plot(gp_pic.X_u, -30, 'rx', 'MarkerSize', 5, 'LineWidth', 2)
plot(x,Eft_pic+2.*sqrt(Varyt_pic),'g--', 'LineWidth', 2)
for i = 1:length(edges)
    plot([edges(i) edges(i)],[-30 35], 'k:')
end
axis tight
caption2 = sprintf('PIC:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f \n s^2_{noise} = %.2f', gp_pic.cf{1}.lengthScale, gp_pic.cf{1}.magnSigma2, gp_pic.cf{2}.lengthScale, gp_pic.cf{2}.magnSigma2, gp_pic.lik.sigma2);
title(caption2)
legend('Data point', 'predicted mean', '2\sigma error', 'inducing input','Location','Northwest')

%========================================================
% PART 4 data analysis with CS+FIC model
%========================================================

% Here we conduct the same analysis as in part 1, but this time we 
% use CS+FIC approximation

% Create the CS+FIC GP structure
if ~exist('ldlchol')
  error('GPstuff:SuiteSparseMissing',...
        ['SuiteSparse is not properly installed. (in BECS try ''use suitesparse'')\n' ...
         'Can not use CS+FIC without SuiteSparse']);
end
gp_csfic = gp_set('type','CS+FIC', 'lik', lik, 'cf', {gpcf1, gpcf2}, 'jitterSigma2', 1e-9, 'X_u', Xu)

% -----------------------------
% --- Conduct the inference ---

% --- MAP estimate using modified Newton algorithm ---

% Now you can choose, if you want to optimize only parameters or
% optimize simultaneously parameters and inducing inputs. Note that
% the inducing inputs are not transformed through logarithm when
% packed

% optimize parameters and inducing inputs
%gp_csfic = gp_set(gp_csfic, 'infer_params', 'covariance+likelihood+inducing');  
% optimize only parameters (default)
%gp_csfic = gp_set(gp_csfic, 'infer_params', 'covariance+likelihood');           

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the scaled conjugate gradient method
gp_csfic=gp_optim(gp_csfic,x,y,'opt',opt);

% Make the prediction
[Eft_csfic, Varft_csfic, lpyt_csfic, Eyt_csfic, Varyt_csfic] = gp_pred(gp_csfic, x, y, x, 'yt', y);

% Plot the solution of FIC
figure
%subplot(4,1,1)
hold on
plot(x,y,'.', 'MarkerSize',7)
plot(x,Eft_csfic,'k', 'LineWidth', 2)
plot(x,Eft_csfic-2.*sqrt(Varyt_csfic),'g--', 'LineWidth', 1)
plot(gp_csfic.X_u, -30, 'rx', 'MarkerSize', 5, 'LineWidth', 2)
plot(x,Eft_csfic+2.*sqrt(Varyt_csfic),'g--', 'LineWidth', 1)
axis tight
caption2 = sprintf('CS+FIC:  l_1= %.2f, s^2_1 = %.2f, \n l_2= %.2f, s^2_2 = %.2f \n s^2_{noise} = %.2f', gp_csfic.cf{1}.lengthScale, gp_csfic.cf{1}.magnSigma2, gp_csfic.cf{2}.lengthScale, gp_csfic.cf{2}.magnSigma2, gp_csfic.lik.sigma2);
title(caption2)
legend('Data point', 'predicted mean', '2\sigma error', 'inducing input','Location','Northwest')

[Eft, Varft, lpyt, Eyt, Varyt] = gp_pred(gp_csfic, x, y, x, 'yt', y);
[Eft1, Varft1] = gp_pred(gp_csfic, x, y, x, 'predcf', 1);
[Eft2, Varft2] = gp_pred(gp_csfic, x, y, x, 'predcf', 2);

figure
set(gcf,'units','centimeters');
set(gcf,'DefaultAxesPosition',[0.08  0.13   0.84   0.85]);
set(gcf,'DefaultAxesFontSize',16)   %6 8
set(gcf,'DefaultTextFontSize',16)   %6 8
hold on
[AX, H1, H2] = plotyy(x, Eft2, x, Eft1+avgy);
set(H2,'LineStyle','--')
set(H2, 'LineWidth', 3)
set(H1,'LineStyle','-')
set(H1, 'LineWidth', 1)

set(AX(2), 'XLim', [-1 559])
set(AX(1), 'XLim', [-1 559])
set(AX(2), 'YLim', [310 380])
set(AX(1), 'YLim', [-5 5])
set(AX(2), 'XTick' ,[0 276 557])
set(AX(2), 'XTicklabel' ,[1958 1981 2004])
set(AX(1), 'XTick' ,[0 276 557])
set(AX(1), 'XTicklabel' ,[1958 1981 2004])
set(AX(2),'YTick',[310 350 380])
set(AX(2),'YTicklabel',[310 350 380])
set(AX(1),'YTick',[-5 0 5])
set(AX(1),'YTicklabel',[-5 0 5])
%set(get(AX(2),'Ylabel'),'String','ppmv')
%set(get(AX(1),'Ylabel'),'String','ppmv') 
set(get(AX(2),'Xlabel'),'String','year')
set(get(AX(1),'Xlabel'),'String','year') 

set(gcf,'pos',[5    3   18  10.7])
set(gcf,'paperunits',get(gcf,'units'))
set(gcf,'paperpos',get(gcf,'pos'))
