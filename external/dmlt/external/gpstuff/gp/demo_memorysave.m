%DEMO_MEMORYSAVE  Demonstration of memory save option in GPstuff
%
%  Description
%    This demo consists of various combinations of covariance functions,
%    likelihoods and full/sparse approximations. This demo is intended for
%    showing how to use memory save option in GPstuff. Unfortunately MATLAB
%    doesn't have memory usage monitoring, so the purpose of the demo is only
%    to test that results are the same whether you use memory save or not
%    and to show that the running times with memory saving are little bit longer
%    (more overhead). Memory saving can be enabled with the following command
%
%       gp = gp_set(..., 'savememory', 'on');
%
%  See also 
%    GP_SET
%

% Copyright (c) 2012 Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Classification data with full GP's

S = which('demo_classific');
L = strrep(S,'demo_classific.m','demodata/synth.tr');
x=load(L);
y=x(:,end);
y = 2.*y-1;
x(:,end)=[];

lik = lik_logit();

pl = prior_t();
pm = prior_sqrtunif();

gpcf = gpcf_neuralnetwork('weightSigma2', [0.9 0.9], 'biasSigma2', 10);
gpcf = gpcf_neuralnetwork(gpcf, 'weightSigma2_prior', pl,'biasSigma2_prior', pm); %

gpcf2 = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 10);
gpcf2 = gpcf_sexp(gpcf2, 'lengthScale_prior', pl,'magnSigma2_prior', pm); %

% GP without memory saving
gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4, 'latent_method', 'EP');
% GP with memory saving option enabled
gp2 = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4, 'latent_method', 'EP', 'savememory', 'on');

fprintf('Classification (EP), Neuralnetwork covariance function without and with memory saving (optimization and prediction)\n');
opt=optimset('TolX',1e-3,'TolFun',1e-3, 'Display', 'off');

% Optimization and prediction without memory saving
tic,gp=gp_optim(gp,x,y,'opt',opt);
[Eft,Varft, lpyt]=gp_pred(gp,x,y,x, 'yt', y);toc

% Optimization and prediction with memory saving
tic,gp2=gp_optim(gp2,x,y,'opt',opt);
[Eft2,Varft2, lpyt2]=gp_pred(gp2,x,y,x, 'yt', y);toc

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Eft==Eft2)); assert(all(Varft==Varft2)); assert(all(lpyt==lpyt2)); 

gp = gp_set('lik', lik, 'cf', gpcf2, 'jitterSigma2', 1e-6, 'latent_method', 'Laplace');
gp2 = gp_set('lik', lik, 'cf', gpcf2, 'jitterSigma2', 1e-6, 'latent_method', 'Laplace', 'savememory', 'on');

fprintf('Classification (Laplace), Squared-Exponential covariance function without and with memory saving (optimization and prediction)\n');

% Optimization and prediction without memory saving
tic,gp=gp_optim(gp,x,y,'opt',opt);
[Eft,Varft, lpyt]=gp_pred(gp,x,y,x, 'yt', y);toc

% Optimization and prediction with memory saving
tic,gp2=gp_optim(gp2,x,y,'opt',opt);
[Eft2,Varft2, lpyt2]=gp_pred(gp2,x,y,x, 'yt', y);toc

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Eft==Eft2)); assert(all(Varft==Varft2)); assert(all(lpyt==lpyt2)); 

% Regression data with sparse approximations

prevstream=setrandstream();
xx=linspace(1,10,901);
x1=logspace(0,1,100);
x1=round(x1*100)-99;
x=xx(x1)';
y=2*sin(4*x)+0.2*randn(size(x));
xt=[1:0.01:14]';
[n,nin] = size(x);

pl = prior_t('s2', 1);
pm = prior_logunif();
pn = prior_logunif();
gpcfse = gpcf_matern32('lengthScale',0.5,'magnSigma2',1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
lik = lik_gaussian('sigma2', 0.1, 'sigma2_prior', pn);
gp = gp_set('lik', lik, 'cf', gpcfse, 'jitterSigma2', 1e-6);
opt=optimset('TolFun',1e-4,'TolX',1e-4);

fprintf('Regression, FIC GP, Matern-3/2 covariance function, without and with memory saving (optimization and prediction)\n')
Xu=round(10+90*rand(18,1))/10; % Random placement

gp_fic = gp_set(gp, 'type','FIC','X_u',Xu,'infer_params','covariance+likelihood+inducing');
gp_fic2 = gp_set(gp, 'type','FIC','X_u',Xu,'infer_params','covariance+likelihood+inducing', 'savememory', 'on');

opt=optimset('TolFun',1e-4,'TolX',1e-4, 'Display', 'off');

% Optimization and prediction without memory saving
tic,gp_fic=gp_optim(gp_fic,x,y,'opt',opt);
[Eft_fic, Varft_fic] = gp_pred(gp_fic, x, y, xt);
Varft_fic = Varft_fic + gp_fic.lik.sigma2;toc

% Optimization and prediction with memory saving
tic,gp_fic2=gp_optim(gp_fic2,x,y,'opt',opt);
[Eft_fic2, Varft_fic2] = gp_pred(gp_fic2, x, y, xt);
Varft_fic2 = Varft_fic2 + gp_fic.lik.sigma2;toc
setrandstream([],prevstream);

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Eft_fic==Eft_fic2)); assert(all(Varft_fic==Varft_fic2));

gpcfse = gpcf_sexp('lengthScale',0.5,'magnSigma2',1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gp = gp_set('lik', lik, 'cf', gpcfse, 'jitterSigma2', 1e-6);

gp_fic = gp_set(gp, 'type','VAR','X_u',Xu,'infer_params','covariance+likelihood+inducing');
gp_fic2 = gp_set(gp, 'type','VAR','X_u',Xu,'infer_params','covariance+likelihood+inducing', 'savememory', 'on');

fprintf('Regression, VAR GP, Squared-Exponential covariance function, without and with memory saving (optimization and prediction)\n')

% Optimization and prediction without memory saving
tic,gp_fic=gp_optim(gp_fic,x,y,'opt',opt);
[Eft_fic, Varft_fic] = gp_pred(gp_fic, x, y, xt);
Varft_fic = Varft_fic + gp_fic.lik.sigma2;toc

% Optimization and prediction with memory saving
tic,gp_fic2=gp_optim(gp_fic2,x,y,'opt',opt);
[Eft_fic2, Varft_fic2] = gp_pred(gp_fic2, x, y, xt);
Varft_fic2 = Varft_fic2 + gp_fic.lik.sigma2;toc
setrandstream([],prevstream);

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Eft_fic==Eft_fic2)); assert(all(Varft_fic==Varft_fic2));


% Spatial data with sparse approximations

S = which('demo_spatial1');
data = load(strrep(S,'demo_spatial1.m','demodata/spatial1.txt'));
x = data(1:200,1:2);
ye = data(1:200,3);
y = data(1:200,4);
dims = [1    60     1    35];
[trindex, Xu] = set_PIC(x, dims, 5, 'corners', 0);
[n,nin] = size(x);

pl = prior_t('s2',10);
pm = prior_sqrtunif();

gpcf1 = gpcf_ppcs3('nin',nin,'lengthScale', 1, 'magnSigma2', 0.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 0.05, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf3 = gpcf_neuralnetwork('weightSigma2', [1 1], 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);

lik = lik_negbin();
% GP without memory save
gp = gp_set('type', 'FIC', 'lik', lik, 'cf', gpcf1, 'X_u', Xu, ...
            'jitterSigma2', 1e-4, 'infer_params', 'covariance+inducing');
% GP with memory saving option enabled
gp2 = gp_set('type', 'FIC', 'lik', lik, 'cf', gpcf1, 'X_u', Xu, ...
            'jitterSigma2', 1e-4, 'infer_params', 'covariance+inducing', 'savememory', 'on');
opt=optimset('TolFun',1e-2,'TolX',1e-2, 'Display', 'off');


fprintf('Spatial process (Laplace), FIC GP, PPCS3 covariance function, without and with memory saving (optimization and prediction)\n');

% Optimization and prediction without memory saving
tic,gp=gp_optim(gp,x,y,'z',ye,'opt',opt);
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye, 'tstind', [1:n]); toc

% Optimization and prediction with memory saving
tic,gp2=gp_optim(gp2,x,y,'z',ye,'opt',opt);
[Ef2, Varf2] = gp_pred(gp2, x, y, x, 'z', ye, 'tstind', [1:n]); toc

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Ef==Ef2)); assert(all(Varf==Varf2));

gp = gp_set('lik', lik, 'cf', gpcf2, 'jitterSigma2', 1e-4, 'latent_method', 'EP');
gp2 = gp_set('lik', lik, 'cf', gpcf2, 'jitterSigma2', 1e-4, 'latent_method', 'EP', 'savememory', 'on');

fprintf('Spatial process (EP), FIC GP, Squared-Exponential covariance function, without and with memory saving (optimization and prediction)\n');

% Optimization and prediction without memory saving
tic,gp=gp_optim(gp,x,y,'z',ye,'opt',opt);
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye, 'tstind', [1:n]); toc

% Optimization and prediction with memory saving
tic,gp2=gp_optim(gp2,x,y,'z',ye,'opt',opt);
[Ef2, Varf2] = gp_pred(gp2, x, y, x, 'z', ye, 'tstind', [1:n]); toc

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Ef==Ef2)); assert(all(Varf==Varf2));

gp = gp_set('lik', lik, 'cf', gpcf3, 'jitterSigma2', 1e-4);
gp2 = gp_set('lik', lik, 'cf', gpcf3, 'jitterSigma2', 1e-4, 'savememory', 'on');

fprintf('Spatial process (Laplace), FIC GP, Neuralnetwork covariance function, without and with memory saving (optimization and prediction)\n');
tic,gp=gp_optim(gp,x,y,'z',ye,'opt',opt);
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye, 'tstind', [1:n]); toc

tic,gp2=gp_optim(gp2,x,y,'z',ye,'opt',opt);
[Ef2, Varf2] = gp_pred(gp2, x, y, x, 'z', ye, 'tstind', [1:n]); toc

% Check that predictions (and optimization) with and without memory saving are same
assert(all(Ef==Ef2)); assert(all(Varf==Varf2));