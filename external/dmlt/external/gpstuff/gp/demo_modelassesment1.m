%DEMO_MODELASSESMENT1 Demonstration for model assessment with WAIC,
% DIC, number of effective parameters and ten-fold cross validation
%                       
%  Description
%    We will consider the regression problem in demo_regression1. 
%    The analysis is conducted with full Gaussian process, and FIC
%    and PIC sparse approximations. The performance of these models
%    are compared by evaluating ten-fold cross validation,
%    leave-one-out cross-validation, WAIC, DIC and the effective
%    number of parameters. The inference will be conducted using
%    maximum a posterior (MAP) estimate for the parameters, via
%    full Markov chain Monte Carlo (MCMC) and with an integration
%    approximation (IA) for the parameters.
%
%    This demo is organised in three parts:
%     1) data analysis with full GP model
%     2) data analysis with FIC approximation
%     3) data analysis with PIC approximation
%
%  See also DEMO_REGRESSION1, DEMO_REGRESSION_SPARSE1

% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010-2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


%========================================================
% PART 1 data analysis with full GP model
%========================================================
disp('Full GP model with Gaussian noise model')

% Load the data
S = which('demo_regression1');
L = strrep(S,'demo_regression1.m','demodata/dat.1');
data=load(L);
x = [data(:,1) data(:,2)];
y = data(:,3);
[n, nin] = size(x);

DIC=repmat(NaN,1,9);DIC2=repmat(NaN,1,9);DIC_latent=repmat(NaN,1,9);
p_eff=repmat(NaN,1,9);p_eff2=repmat(NaN,1,9);p_eff_latent=repmat(NaN,1,9);p_eff_latent2=repmat(NaN,1,9);

% ---------------------------
% --- Construct the model ---
gpcf = gpcf_sexp('lengthScale', [1 1], 'magnSigma2', 0.2^2);
lik = lik_gaussian('sigma2', 0.2^2);
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);

% -----------------------------
% --- Conduct the inference ---
%
% We will make the inference first by finding a maximum a posterior estimate 
% for the parameters via gradient based optimization. After this we will
% perform an extensive Markov chain Monte Carlo sampling for the parameters.
% 
disp('MAP estimate for the parameters')

% --- MAP estimate ---
%     (see gp_optim for more details)

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters, DIC and WAIC with focus on
% latent variables. 
models{1} = 'full_MAP';
p_eff_latent = gp_peff(gp, x, y);
% For easier comparison to other methods, compute mean log
% predictive density (mlpd) instead of deviance (-2n*mlpd)
[DIC_latent, p_eff_latent2] = gp_dic(gp, x, y, 'focus', 'latent', 'output', 'mlpd');
WAICV(1) = gp_waic(gp,x,y);
WAICG(1) = gp_waic(gp,x,y, 'method', 'G');


% Evaluate the 10-fold cross-validation results.
disp('MAP estimate for the parameters - k-fold-CV')
cvres =  gp_kfcv(gp, x, y, 'display', 'fold');
mlpd_cv(1) = cvres.mlpd_cv;
rmse_cv(1) = cvres.mrmse_cv;

disp('MAP estimate for the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gp, x, y);
mlpd_loo(1) = mean(lpy);
rmse_loo(1) = sqrt(mean((y-Ef).^2));

% --- MCMC approach ---
disp('MCMC integration over the parameters')

% Do the sampling (this takes about 1 minute)
[rfull,g,opt] = gp_mc(gp, x, y, 'nsamples', 220, 'display', 20);

% After sampling delete the burn-in and thin the sample chain
rfull = thin(rfull, 21, 2);

% Evaluate the effective number of parameters, DIC and WAIC. 
models{2} = 'full_MCMC';
% For easier comparison to other methods, compute mean log
% predictive density (mlpd) instead of deviance (-2n*mlpd)
[DIC(2), p_eff(2)] =  gp_dic(rfull, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(2), p_eff2(2)] =  gp_dic(rfull, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(2) = gp_waic(rfull,x,y);
WAICG(2) = gp_waic(rfull,x,y, 'method', 'G');

% Evaluate the 10-fold cross validation results. 
%
% We reduce the number of samples so that the sampling takes less time. 
% 50 is too small sample size, though, and for reliable results the 10-CV 
% should be run with larger sample size.
disp('MCMC integration over the parameters - k-fold-CV')
opt.nsamples= 50; 
cvres =  gp_kfcv(gp, x, y, 'inf_method', 'MCMC', 'opt', opt, 'rstream', 1, 'display', 'fold');
mlpd_cv(2) = cvres.mlpd_cv;
rmse_cv(2) = cvres.mrmse_cv;

disp('MCMC integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(rfull, x, y);
mlpd_loo(2) = mean(lpy);
rmse_loo(2) = sqrt(mean((y-Ef).^2));

% --- Integration approximation approach ---
disp('Grid integration over the parameters')

gp_array = gp_ia(gp, x, y, 'int_method', 'grid');

models{3} = 'full_IA'; 
% For easier comparison to other methods, compute mean log
% predictive density (mlpd) instead of deviance (-2n*mlpd)
[DIC(3), p_eff(3)] =  gp_dic(gp_array, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(3), p_eff2(3)] =  gp_dic(gp_array, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(3) = gp_waic(gp_array,x,y);
WAICG(3) = gp_waic(gp_array,x,y, 'method', 'G');

% Then the 10 fold cross-validation.
disp('Grid integration over the parameters - k-fold-CV')
clear opt
opt.int_method = 'grid';
cvres = gp_kfcv(gp, x, y, 'inf_method', 'IA', 'opt', opt, 'display', 'fold');
mlpd_cv(3) = cvres.mlpd_cv;
rmse_cv(3) = cvres.mrmse_cv;

disp('Grid integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gp_array, x, y);
mlpd_loo(3) = mean(lpy);
rmse_loo(3) = sqrt(mean((y-Ef).^2));

%========================================================
% PART 2 data analysis with FIC GP
%========================================================
disp('GP with FIC sparse approximation')

% ---------------------------
% --- Construct the model ---

% Here we conduct the same analysis as in part 1, but this time we 
% use FIC approximation

% Initialize the inducing inputs in a regular grid over the input space
[u1,u2]=meshgrid(linspace(-1.8,1.8,6),linspace(-1.8,1.8,6));
X_u = [u1(:) u2(:)];

% Create the FIC GP structure
gp_fic = gp_set('type', 'FIC', 'lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-6, 'X_u', X_u)

% -----------------------------
% --- Conduct the inference ---

% --- MAP estimate using scaled conjugate gradient algorithm ---
disp('MAP estimate for the parameters')

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp_fic=gp_optim(gp_fic,x,y,'opt',opt,'optimf',@fminlbfgs);

% Evaluate the effective number of parameters and DIC with focus on
% latent variables.
models{4} = 'FIC_MAP';
p_eff_latent(4) = gp_peff(gp_fic, x, y);
[DIC_latent(4), p_eff_latent2(4)] = gp_dic(gp_fic, x, y, 'focus', 'latent', 'output', 'mlpd');
WAICV(4) = gp_waic(gp_fic,x,y);
WAICG(4) = gp_waic(gp_fic,x,y, 'method', 'G');

% Evaluate the 10-fold cross validation results. 
disp('MAP estimate for the parameters - k-fold-CV')
cvres = gp_kfcv(gp_fic, x, y, 'display', 'fold');
mlpd_cv(4) = cvres.mlpd_cv;
rmse_cv(4) = cvres.mrmse_cv;

disp('MAP estimate for the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gp_fic, x, y);
mlpd_loo(4) = mean(lpy);
rmse_loo(4) = sqrt(mean((y-Ef).^2));

% --- MCMC approach ---
% (the inducing inputs are fixed)
disp('MCMC integration over the parameters')

% Do the sampling (this takes about 1 minute)
rfic = gp_mc(gp_fic, x, y, 'nsamples', 220, 'display', 20);

% After sampling we delete the burn-in and thin the sample chain
rfic = thin(rfic, 21, 2);

% Evaluate the effective number of parameters, DIC and WAIC. Note that 
% the effective number of parameters as a second output, but here 
% we use explicitly the gp_peff function
models{5} = 'FIC_MCMC'; 
[DIC(5), p_eff(5)] =  gp_dic(rfic, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(5), p_eff2(5)] =  gp_dic(rfic, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(5) = gp_waic(rfic,x,y);
WAICG(5) = gp_waic(rfic,x,y, 'method', 'G');

% We reduce the number of samples so that the sampling takes less time. 
% 50 is too small sample size, though, and for reliable results the 10-CV 
% should be run with larger sample size. We also set the save option to 0.
clear opt
opt.nsamples= 50; opt.display=20; 
disp('MCMC integration over the parameters - k-fold-CV')
cvres = gp_kfcv(gp_fic, x, y, 'inf_method', 'MCMC', 'opt', opt, 'display', 'fold');
mlpd_cv(5) = cvres.mlpd_cv;
rmse_cv(5) = cvres.mrmse_cv;

disp('MCMC integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(rfic, x, y);
mlpd_loo(5) = mean(lpy);
rmse_loo(5) = sqrt(mean((y-Ef).^2));

% --- Integration approximation approach ---
disp('Grid integration over the parameters')
gpfic_array = gp_ia(gp_fic, x, y, 'int_method', 'grid');

models{6} = 'FIC_IA'; 
[DIC(6), p_eff(6)] =  gp_dic(gpfic_array, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(6), p_eff2(6)] =  gp_dic(gpfic_array, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(6) = gp_waic(gpfic_array,x,y);
WAICG(6) = gp_waic(gpfic_array,x,y, 'method', 'G');

% Then the 10 fold cross-validation.
disp('Grid integration over the parameters - k-fold-CV')
clear opt
opt.int_method = 'grid';
cvres = gp_kfcv(gp_fic, x, y, 'inf_method', 'IA', 'opt', opt, 'display', 'fold');
mlpd_cv(6) = cvres.mlpd_cv;
rmse_cv(6) = cvres.mrmse_cv;

disp('Grid integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gpfic_array, x, y);
mlpd_loo(6) = mean(lpy);
rmse_loo(6) = sqrt(mean((y-Ef).^2));

%========================================================
% PART 3 data analysis with PIC approximation
%========================================================
disp('GP with PIC sparse approximation')

[u1,u2]=meshgrid(linspace(-1.8,1.8,6),linspace(-1.8,1.8,6));
X_u = [u1(:) u2(:)];

% Initialize test points
[p1,p2]=meshgrid(-1.8:0.1:1.8,-1.8:0.1:1.8);
p=[p1(:) p2(:)];

% set the data points into clusters. Here we construct two cell arrays. 
%  trindex  contains the block index vectors for training data. That is 
%           x(trindex{i},:) and y(trindex{i},:) belong to the i'th block.
b1 = [-1.7 -0.8 0.1 1 1.9];
mask = zeros(size(x,1),size(x,1));
trindex={}; 
for i1=1:4
  for i2=1:4
    ind = 1:size(x,1);
    ind = ind(: , b1(i1)<=x(ind',1) & x(ind',1) < b1(i1+1));
    ind = ind(: , b1(i2)<=x(ind',2) & x(ind',2) < b1(i2+1));
    trindex{4*(i1-1)+i2} = ind';
  end
end

% Create the PIC GP structure and set the inducing inputs and block indexes
gpcf = gpcf_sexp('lengthScale', [1 1], 'magnSigma2', 0.2^2);
lik = lik_gaussian('sigma2', 0.2^2);

gp_pic = gp_set('type', 'PIC', 'lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-6, 'X_u', X_u);
gp_pic = gp_set(gp_pic, 'tr_index', trindex);

% -----------------------------
% --- Conduct the inference ---

% --- MAP estimate using scaled conjugate gradient algorithm ---
disp('MAP estimate for the parameters')

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the BFGS quasi-Newton method
gp_pic=gp_optim(gp_pic,x,y,'opt',opt,'optimf',@fminlbfgs);

models{7} = 'PIC_MAP';
p_eff_latent(7) = gp_peff(gp_pic, x, y);
[DIC_latent(7), p_eff_latent2(7)] = gp_dic(gp_pic, x, y, 'focus', 'latent', 'output', 'mlpd');
WAICV(7) = gp_waic(gp_pic, x, y);
WAICG(7) = gp_waic(gp_pic, x, y, 'method', 'G');

% Evaluate the 10-fold cross validation results. 
disp('MAP estimate for the parameters - k-fold-CV')
cvres = gp_kfcv(gp_pic, x, y, 'display', 'fold');
mlpd_cv(7) = cvres.mlpd_cv;
rmse_cv(7) = cvres.mrmse_cv;

disp('MAP estimate for the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gp_pic, x, y);
mlpd_loo(7) = mean(lpy);
rmse_loo(7) = sqrt(mean((y-Ef).^2));

% --- MCMC approach ---
disp('MCMC integration over the parameters')

% Do the sampling (this takes about 1 minute)
rpic = gp_mc(gp_pic, x, y, 'nsamples', 220, 'display', 20);

% After sampling we delete the burn-in and thin the sample chain
rpic = rmfield(rpic, 'tr_index');
rpic = thin(rpic, 21, 2);
rpic.tr_index = trindex;

% Evaluate the effective number of parameters and DIC. Note that 
% the effective number of parameters as a second output, but here 
% we use explicitly the gp_peff function
models{8} = 'PIC_MCMC'; 
[DIC(8), p_eff(8)] =  gp_dic(rpic, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(8), p_eff2(8)] =  gp_dic(rpic, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(8) = gp_waic(rpic, x, y);
WAICG(8) = gp_waic(rpic, x, y, 'method', 'G');

% We reduce the number of samples so that the sampling takes less time. 
% 50 is too small sample size, though, and for reliable results the 10-CV 
% should be run with larger sample size. We also set the save option to 0.
clear opt
opt.nsamples= 50; opt.display=20;
disp('MCMC integration over the parameters - k-fold-CV')
cvres = gp_kfcv(gp_pic, x, y, 'inf_method', 'MCMC', 'opt', opt, 'display', 'fold');
mlpd_cv(8) = cvres.mlpd_cv;
rmse_cv(8) = cvres.mrmse_cv;

disp('MCMC integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(rpic, x, y);
mlpd_loo(8) = mean(lpy);
rmse_loo(8) = sqrt(mean((y-Ef).^2));

% --- Integration approximation approach ---
disp('Grid integration over the parameters')

gppic_array = gp_ia(gp_pic, x, y, 'int_method', 'grid');

models{9} = 'PIC_IA'; 
[DIC(9), p_eff(9)] =  gp_dic(gppic_array, x, y, 'focus', 'param', 'output', 'mlpd');
[DIC2(9), p_eff2(9)] =  gp_dic(gppic_array, x, y, 'focus', 'all', 'output', 'mlpd');
WAICV(9) = gp_waic(gppic_array, x, y);
WAICG(9) = gp_waic(gppic_array, x, y, 'method', 'G');

% Then the 10 fold cross-validation.
disp('Grid integration over the parameters - k-fold-CV')
clear opt
opt.int_method = 'grid';
cvres = gp_kfcv(gp_pic, x, y, 'inf_method', 'IA', 'opt', opt, 'display', 'fold');
mlpd_cv(9) = cvres.mlpd_cv;
rmse_cv(9) = cvres.mrmse_cv;

disp('Grid integration over the parameters - LOO-CV')
[Ef,Varf,lpy] =  gp_loopred(gppic_array, x, y);
mlpd_loo(9) = mean(lpy);
rmse_loo(9) = sqrt(mean((y-Ef).^2));

%========================================================
% PART 4 Print the results
%========================================================
disp('Summary of the results')
S = '      ';
for i = 1:length(models)
    S = [S '  ' models{i}];
end

S = sprintf([S '\n CV-mlpd  %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], mlpd_cv);
S = sprintf([S '\n CV-rmse  %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], rmse_cv);
S = sprintf([S '\n LOO-mlpd %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], mlpd_loo);
S = sprintf([S '\n LOO-rmse %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], rmse_loo);
S = sprintf([S '\n ']);
S = sprintf([S '\n WAIC_V   %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], WAICV);
S = sprintf([S '\n WAIC_G   %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], WAICG);
S = sprintf([S '\n ']);
S = sprintf([S '\n DIC_h    %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], DIC);
S = sprintf([S '\n DIC_a    %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], DIC2);
S = sprintf([S '\n DIC_l    %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], DIC_latent);
S = sprintf([S '\n peff_h   %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], p_eff);
S = sprintf([S '\n peff_a   %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], p_eff2);
S = sprintf([S '\n peff_l   %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], p_eff_latent);
S = sprintf([S '\n peff_l2  %5.2f     %5.2f     %5.2f     %5.2f    %5.2f     %5.2f    %5.2f    %5.2f    %5.2f'], p_eff_latent2);
S = sprintf([S '\n ']);
S = sprintf([S '\n ']);
S = sprintf([S '\n The notation is as follows:']);
S = sprintf([S '\n CV-mlpd  = mean log predictive density from the 10-fold CV. ']);
S = sprintf([S '\n CV-rmse  = root mean squared error from the 10-fold LOO-CV. ']);
S = sprintf([S '\n LOO-mlpd = mean log predictive density from the LOO-CV. ']);
S = sprintf([S '\n LOO-rmse = root mean squared error from the 10-fold CV. ']);
S = sprintf([S '\n WAIC_V   = WAIC via variance method ']);
S = sprintf([S '\n WAIC_G   = WAIC via Gibbs training utility method ']);
S = sprintf([S '\n DIC_h    = DIC with focus on parameters. ']);
S = sprintf([S '\n DIC_a    = DIC with focus on parameters and latent variables (all). ']);
S = sprintf([S '\n DIC_l    = DIC with focus on latent variables. ']);
S = sprintf([S '\n peff_h   = effective number of parameters (latent variables marginalized). ']);
S = sprintf([S '\n peff_a   = effective number of parameters and latent variables. ']);
S = sprintf([S '\n peff_l   = effective number of latent variables evaluated with gp_peff. ']);
S = sprintf([S '\n peff_l2  = effective number of latent variables evaluated with gp_dic. ']);
S = sprintf([S '\n '])
