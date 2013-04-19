%DEMO_INPUTDEPENDENTNOISE Demonstration of input dependent-noise
%                         model using Gaussian process prior
%
% Description
%       Uses toy data sets to demonstrate how inferring
%       heteroscedastic noise with input dependent noise model
%       differs from standard noise models (Gaussian, Student-t).
%

% Copyright (c) Ville Tolvanen 2011-2012
% 
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


%=================================
% 1D Demonstration
%=================================

% Initialize random stream
prevstream=setrandstream(0);

% Create toy data
% x = 100*rand([40 1]);
n = 500;
x=linspace(-100,200,n)';
f1 = [5.*sin(-3+0.2.*x(1:ceil(0.23*n))); 20*sin(0.1*x(ceil(0.23*n)+1:ceil(0.85*n))); 5.*sin(2.8+0.2.*x(ceil(0.85*n)+1:end))];
f2 = 100*norm_pdf(x,110,20) + 100*norm_pdf(x,-10,20);
sigma2 = 0.5;

x=x-mean(x); x=x./std(x);
f1 = f1-mean(f1); f1=f1./std(f1);

y = f1 + sqrt((sigma2.*exp(f2))).*randn(size(x));
yt = f1(1:5:end);
xt = x(1:5:end);
nt = size(xt,1);
x=x(:); y=y(:); xt=xt(:);

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_t('s2',10); 
gpcf1 = gpcf_sexp('lengthScale', 0.5, 'magnSigma2', 0.5, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_exp('lengthScale', 1, 'magnSigma2', 0.1, ...
                 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure. Don't set prior for sigma2 if covariance
% function magnitude for noise process has a prior.
lik = lik_inputdependentnoise('sigma2', 0.1, 'sigma2_prior', prior_fixed());

% NOTE! if multiple covariance functions per latent is used, define
% gp.comp_cf as follows:
% gp = gp_set(..., 'comp_cf' {[1 2] [5 6]};
gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-9, 'comp_cf', {[1] [2]});

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');
% For more complex problems, maxiter in latent_opt should be increased.
% gp.latent_opt.maxiter=1e6;

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Derivativecheck','off');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

% make prediction to the data points
[Ef, Varf,lpyt] = gp_pred(gp, x, y, xt, 'yt', yt);
Ef11=Ef(1:nt);Ef12=Ef(nt+1:end);
Varf11=diag(Varf(1:nt,1:nt));
%prctmus = gp_predprctmu(gp, x, y, xt);
prctmus=[Ef11-1.645*sqrt(Varf11) Ef11 Ef11+1.645*sqrt(Varf11)];
fprintf('mlpd inputdependentnoise: %.2f\n', mean(lpyt));

% Gaussian for comparison
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Derivativecheck','off');
lik2 = lik_gaussian();
gp2 = gp_set('lik', lik2, 'cf', gpcf1, 'jitterSigma2', 1e-9);
gp2 = gp_optim(gp2,x,y,'opt',opt);
[Ef2, Varf2, lpyt2] = gp_pred(gp2, x, y, xt,'yt',yt);
prctmus2 = gp_predprctmu(gp2, x, y, xt);
fprintf('mlpd gaussian: %.2f\n', mean(lpyt2));

% Student-t for comparison
lik=lik_t();
gp3=gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-9);
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Derivativecheck','off');
gp3 = gp_set(gp3, 'latent_method', 'Laplace');
gp3=gp_optim(gp3,x,y,'opt',opt);
[Ef3, Varf3,lpyt3] = gp_pred(gp3, x, y, xt, 'yt', yt);
prctmus3 = gp_predprctmu(gp3, x, y, xt);
fprintf('mlpd student-t: %.2f\n', mean(lpyt3));

figure;
% plot mean and 5% and 95% quantiles
subplot(3,1,1)
plot(xt,Ef11,'b',xt,prctmus(:,1),'r',xt,prctmus(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Input dependent noise model');
legend('Mean','5%','95%','True',2)

% Compare to Gaussian with homoscedastic scale
subplot(3,1,2),
plot(xt, Ef2,'b',xt,prctmus2(:,1),'r',xt,prctmus2(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Gaussian noise model')

% Compare to Student-t with homoscedastic scale
subplot(3,1,3)
plot(xt, Ef3,'b',xt,prctmus3(:,1),'r',xt,prctmus3(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Student-t noise model')

figure
s2=gp.lik.sigma2;
plot(xt, s2.*exp(Ef12), '-b',x, sigma2.*exp(f2), '-k', xt, s2.*exp(Ef12 + 1.96.*sqrt(diag(Varf(nt+1:end, nt+1:end)))), '-r', xt,s2.*exp(Ef12 - 1.96.*sqrt(diag(Varf(nt+1:end, nt+1:end)))), '-r')
legend('Predicted noise variance', 'Real noise variance','95% CI',2);

%====================================
% 2D Demonstration
%====================================
setrandstream(0);

% Create data from two 2 dimensional gaussians
nt=10;
n=700;
x=[-3+6*rand(0.25*n,1) -3+6*rand(0.25*n,1);-1.5+3*rand(0.75*n,1) -1.5+3*rand(0.75*n,1)];
mu=[0 0];
S=[0.2 0;0 0.2];
sigma2=0.1;
[x1,x2]=meshgrid(linspace(-1,2,nt), linspace(-1,2,nt));
xt=[x1(:) x2(:)];
f1t=10*mnorm_pdf(xt,mu,S) + 20*mnorm_pdf(xt,mu+[1.5 1.5], [0.5 0;0 0.5]);
f2t=20*mnorm_pdf(xt,mu, [0.9 0;0 0.9]);
yt=f1t;
f1 = 10*mnorm_pdf(x,mu, S) + 20*mnorm_pdf(x,mu+[1.5 1.5], [0.5 0;0 0.5]);
f2 = 20*mnorm_pdf(x,mu, [0.9 0;0 0.9]);

y=f1+randn(size(x,1),1).*sqrt(sigma2.*exp(f2));

pl = prior_logunif();
pm = prior_logunif(); 
gpcf1 = gpcf_sexp('lengthScale', [1 1.01], 'magnSigma2', 1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('lengthScale', [1 1.01], 'magnSigma2', 0.1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);

lik=lik_inputdependentnoise('sigma2', 0.1, 'sigma2_prior', prior_fixed());

gp=gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-6, 'comp_cf', {[1] [2]});
gp = gp_set(gp, 'latent_method', 'Laplace');
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Display','iter','MaxIter',100,'Derivativecheck','off');
gp=gp_optim(gp,x,y,'opt',opt);
% Increase maxiter for predictions in case of slow convergence
gp.latent_opt.maxiter=1e6;
[Ef,Varf,lpyt]=gp_pred(gp,x,y,xt, 'yt',yt);
fprintf('mlpd inputdependentnoise: %.2f\n', mean(lpyt));

% Gaussian for comparison
lik2 = lik_gaussian('sigma2', sigma2);
gp2 = gp_set('lik', lik2, 'cf', gpcf1, 'jitterSigma2', 1e-6);
gp2 = gp_optim(gp2,x,y,'opt',opt);
[Ef2,Varf2,lpyt2]=gp_pred(gp2,x,y,xt,'yt',yt);
fprintf('mlpd gaussian: %.2f\n', mean(lpyt2));

% Student-t for comparison
lik3=lik_t('sigma2', sigma2);
gp3=gp_set('lik', lik3, 'cf', gpcf1, 'jitterSigma2', 1e-6, 'latent_method', 'Laplace');
gp3=gp_optim(gp3,x,y,'opt',opt);
[Ef3,Varf3,lpyt3]=gp_pred(gp3,x,y,xt,'yt',yt);
fprintf('mlpd student-t: %.2f\n', mean(lpyt3));

s2=gp.lik.sigma2;
figure
subplot(3,1,1),mesh(x1,x2,reshape(f1t,size(x1))),hold on, plot3(xt(:,1),xt(:,2), Ef(1:size(xt,1)), '*')
title('Input dependent noise model');
colormap hsv, alpha(.4)
subplot(3,1,2),mesh(x1,x2,reshape(f1t,size(x1))),hold on, plot3(xt(:,1),xt(:,2), Ef2(1:size(xt,1)), '*');
colormap hsv, alpha(.4)
title('Gaussian noise model');
subplot(3,1,3),mesh(x1,x2,reshape(f1t,size(x1))),hold on, plot3(xt(:,1),xt(:,2), Ef3(1:size(xt,1)), '*');
colormap hsv, alpha(.4)
title('Student-t noise model');

figure
mesh(x1,x2,sigma2.*exp(reshape(f2t,size(x1)))),hold on, plot3(xt(:,1),xt(:,2), s2.*exp(Ef(size(xt,1)+1:end)), '*'); 
title('Real noise versus predicted noise');
colormap hsv, alpha(.4)

%============================================
% Demonstration with homoscedastic noise
%============================================
setrandstream(0);

% Create data
n =200;
nt = 200;
x = linspace(-100,200, n)';
xt = linspace(-100,200, n)';
f1 = [5.*sin(-3+0.2.*x(1:ceil(0.23*n))); 20*sin(0.1*x(ceil(0.23*n)+1:ceil(0.85*n))); 5.*sin(2.8+0.2.*x(ceil(0.85*n)+1:end))];
sigma2 = 1;

x=x-mean(x); x=x./std(x);
xt=xt-mean(xt); xt=xt./std(xt);
f1 = f1-mean(f1); f1=f1./std(f1);
f2 = zeros(size(f1));

y = f1 + sqrt(sigma2).*randn(size(x));yt=f1;
x=x(:); y=y(:); xt=xt(:);

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_t('s2',10); 
gpcf1 = gpcf_sexp('lengthScale', 0.5, 'magnSigma2', 0.5, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('lengthScale', 1, 'magnSigma2',0.1, ...
                  'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the the model
lik = lik_inputdependentnoise('sigma2', 0.1, 'sigma2_prior', prior_fixed());
gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-9, 'comp_cf', {[1] [2]});

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');
opt=optimset('TolFun',1e-4,'TolX',1e-4);

% if flat priors are used, there might be need to increase
% gp.latent_opt.maxiter for laplace algorithm to converge properly

% gp.latent_opt.maxiter=1e6;

gp=gp_optim(gp,x,y,'opt',opt);

[Ef, Varf,lpyt] = gp_pred(gp, x, y, xt,'yt',yt);
Ef11=Ef(1:nt);Ef12=Ef(nt+1:end);
Varf11=diag(Varf(1:nt,1:nt));
%prctmus = gp_predprctmu(gp, x, y, xt);
prctmus=[Ef11-1.645*sqrt(Varf11) Ef11 Ef11+1.645*sqrt(Varf11)];
fprintf('mlpd inputdependentnoise: %.2f\n', mean(lpyt));

% Gaussian for comparison
lik2 = lik_gaussian();
gp2 = gp_set('lik', lik2, 'cf', gpcf1, 'jitterSigma2', 1e-6);
gp2 = gp_optim(gp2,x,y,'opt',opt);
[Ef2, Varf2, lpyt2] = gp_pred(gp2, x, y, xt,'yt',yt);
prctmus2 = gp_predprctmu(gp2, x, y, xt);
fprintf('mlpd gaussian: %.2f\n', mean(lpyt2));

% Student-t for comparison
lik=lik_t();
gp3=gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);
gp3 = gp_set(gp3, 'latent_method', 'Laplace');
gp3=gp_optim(gp3,x,y,'opt',opt);
[Ef3, Varf3,lpyt3] = gp_pred(gp3, x, y, xt,'yt',yt);
prctmus3 = gp_predprctmu(gp3, x, y, xt);
fprintf('mlpd student-t: %.2f\n', mean(lpyt3));

figure;
% plot mean and 5% and 95% quantiles
subplot(3,1,1)
plot(xt,Ef11,'b',xt,prctmus(:,1),'r',xt,prctmus(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Input dependent noise model');
legend('Mean','5%','95%','True',2)

% Compare to Gaussian with homoscedastic scale
subplot(3,1,2),
plot(xt, Ef2,'b',xt,prctmus2(:,1),'r',xt,prctmus2(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Gaussian noise model')

% Compare to Student-t with homoscedastic scale
subplot(3,1,3)
plot(xt, Ef3,'b',xt,prctmus3(:,1),'r',xt,prctmus3(:,3),'r', x, f1, 'k')
ylim([-3 3]), title('Student-t noise model')

figure
plot(xt, s2.*exp(Ef12), '-b',x, sigma2.*exp(f2), '-k', xt, s2.*exp(Ef12 + 1.96.*sqrt(diag(Varf(nt+1:end, nt+1:end)))), '-r', xt,s2.*exp(Ef12 - 1.96.*sqrt(diag(Varf(nt+1:end, nt+1:end)))), '-r')
ylim([0 2.5])
legend('Predicted noise variance', 'Real noise variance','95% CI',2);
setrandstream(prevstream);