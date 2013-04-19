%DEMO_HURDLE Demonstration of Logit Negative-binomial hurdle model
%            using Gaussian process prior
%
%  Description
%    Hurdle models can be used to model excess number of zeros
%    compared to usual Poisson and negative binomial count models. 
%    Hurdle models assume a two-stage process, where the first
%    process determines whether the count is larger than zero, and
%    the second process determines the non-zero count. Both
%    processes are given a zero mean Gaussian process prior. The
%    two stage model formulation makes it possible to make
%    inference for the two latent processes separately.
%
%    In this demo we construct logit negative binomial hurdle model
%    by using logit and zero truncated negative binomial models. 
%    The posterior inference is made with parallel-EP
%    approximation.
%
%    See also  DEMO_SPATIAL2, DEMO_CLASSIFIC1

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari
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

% for zero-process
yz=double(y>0)*2-1;
% for count-process
ci=find(y>0);
yc=y(ci);
xc=x(ci,:);
yec=ye(ci,:);

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_sqrtunif();
cf = gpcf_matern32('lengthScale', 5, 'magnSigma2', 0.05, ...
                   'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the zero part
likz=lik_probit();
gpz=gp_set('lik',likz,'cf',cf,'jitterSigma2',1e-6,'latent_method','EP','latent_opt',struct('parallel','on'));
% Create the count part
likc=lik_negbinztr();
gpc=gp_set('lik',likc,'cf',cf,'jitterSigma2',1e-6,'latent_method','EP','latent_opt',struct('parallel','on'));
% Note that although both parts use 'cf', they have separate hyperparameters

% Set the options for the quasi-Newton optimization
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gpz=gp_optim(gpz,x,yz,'opt',opt,'optimf',@fminlbfgs);
gpc=gp_optim(gpc,xc,yc,'z',yec,'opt',opt,'optimf',@fminlbfgs);

% make prediction to the data points
[Efz, Varfz] = gp_pred(gpz, x, yz, x);
[Efc, Varfc] = gp_pred(gpc, xc, yc, xc, 'z', yec);

% Define help parameters for plotting
xii=sub2ind([60 35],x0(:,2),x0(:,1));
[X1,X2]=meshgrid(1:35,1:60);

% Plot the figures
figure
subplot(1,2,1)
G=repmat(NaN,size(X1));
G(xii)=Efz(:);
pcolor(X1,X2,G),shading flat
colorbar
axis equal
axis([0 35 0 60])
title('Posterior mean of latent zero process')

subplot(1,2,2)
G=repmat(NaN,size(X1));
G(xii(ci))=Efc(:);
pcolor(X1,X2,G),shading flat
colorbar
axis equal
axis([0 35 0 60])
title('Posterior mean of latent count process')
