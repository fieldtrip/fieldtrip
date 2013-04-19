%DEMO_SURVIVAL_COXPH  Survival model using Cox proportional model 
%
%  Description 
%    Survival model using Cox proportional model with a piecewise
%    log-constant baseline hazard. The hazard rate is 
%   
%       h(t) = h_0(t)*exp(f),
%
%    where the baseline hazard is assumed to piecewise log-constant. 
%
%    The inference is conducted via Laplace, where we find
%    Gaussian approximation for p(f| th, data), where th is the
%    maximum a posterior (MAP) estimate for the parameters.
%    
%    The censoring indicator ye is
%    
%      ye = 0 for uncensored event
%      ye = 1 for right censored event.
%
%    If survival times y for n observation are given as nx2 matrix with
%    entry times into follow-up in the first column and exit times from
%    follow-up in the second column, left truncated right censored
%    modelling is possible, for instance, in cases where age is wanted to
%    be set as a baseline hazard.  
%
%    Example data set is leukemia survival data in Northwest England
%    presented in (Henderson, R., Shimakura, S., and Gorst, D. (2002).
%    Modeling spatial variation in leukemia survival data. Journal of the
%    American Statistical Association, 97:965–972). Data set was downloaded
%    from http://www.math.ntnu.no/%7Ehrue/r-inla.org/examples/leukemia/leuk.dat
%
%  See also  DEMO_SURVIVAL_WEIBULL

% Copyright (c) 2011 Jaakko Riihimäki

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% First load data
S = which('demo_survival_weibull');
L = strrep(S,'demo_survival_weibull.m','demodata/leukemia.txt');
leukemiadata=load(L);

% leukemiadata consists of:
% 'time', 'cens', 'xcoord', 'ycoord', 'age', 'sex', 'wbc', 'tpi', 'district'

% survival times
y=leukemiadata(:,1);
% scale survival times
y=y/max(y);

ye=1-leukemiadata(:,2); % event indicator, ye = 0 for uncensored event
                  %                        ye = 1 for right censored event

% choose (for example) 'age', 'sex', 'wbc', and 'tpi' covariates
x0=leukemiadata(:,5:8);
x=x0;
% normalize continuous covariates 
x(:,[1 3:4])=bsxfun(@rdivide,bsxfun(@minus,x0(:,[1 3:4]),mean(x0(:,[1 3:4]),1)),std(x0(:,[1 3:4]),1));

[n, nin]=size(x);

% number of time intervals
ntime=50;
% create finite partition of time axis
S=linspace(0,max(y)+0.001,ntime+1);

% Create the covariance functions
pl = prior_t('s2',1, 'nu', 4);
pm = prior_t('s2',1, 'nu', 4); 

% covariance for hazard function
gpcfh = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% covariance for proportional part
gpcf = gpcf_sexp('lengthScale', ones(1,size(x,2)), 'magnSigma2', 1.2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
lik = lik_coxph('S', S);

% NOTE! if multiple covariance functions per latent is used, define
% gp.comp_cf as follows:
% gp = gp_set(..., 'comp_cf' {[1 2] [5 6]};
% where [1 2] are for hazard function, and [5 6] for proportional part
gp = gp_set('lik', lik, 'cf', {gpcfh gpcf}, 'jitterSigma2', 1e-6, 'comp_cf', {1 2});

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');

opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gp=gp_optim(gp,x,y,'z',ye,'opt',opt);

% Make prediction
xt1=zeros(200,nin); xt1(:,2)=1;
xt2=zeros(200,nin); xt2(:,2)=-1;

xt1(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt2(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt01(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);
xt02(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);


[Ef1, Covf1] = gp_pred(gp, x, y, xt1, 'z', ye);
[Ef2, Covf2] = gp_pred(gp, x, y, xt2, 'z', ye);
Varf1 = diag(Covf1);
Varf2 = diag(Covf2);


figure, hold on, set(gcf, 'color', 'w'),
plot(gp.lik.xtime, Ef1(1:ntime),'k', 'linewidth', 3)
plot(gp.lik.xtime, Ef1(1:ntime)+1.96*sqrt(Varf1(1:ntime)),'--k', 'linewidth', 2)
plot(gp.lik.xtime, Ef1(1:ntime)-1.96*sqrt(Varf1(1:ntime)),'--k', 'linewidth', 2)
title('log-baseline hazard (follow-up time)')
xlabel('time')

% Compute posterior mean and 95% credible intervals of latent function as a
% function of age and sex
col1=ones(1,3)*0.7;
col2=ones(1,3)*0.3;
figure, hold on, set(gcf, 'color', 'w'),
plot(xt01(:,1), Ef1(ntime+1:end), 'color', col1, 'linewidth', 3)
plot(xt01(:,1), Ef1(ntime+1:end)+1.96*sqrt(Varf1(ntime+1:end)), '--', 'color', col1, 'linewidth', 2)
plot(xt01(:,1), Ef1(ntime+1:end)-1.96*sqrt(Varf1(ntime+1:end)), '--', 'color', col1, 'linewidth', 2)

plot(xt02(:,1), Ef2(ntime+1:end), 'color', col2, 'linewidth', 3)
plot(xt02(:,1), Ef2(ntime+1:end)+1.96*sqrt(Varf2(ntime+1:end)), '--', 'color', col2, 'linewidth', 2)
plot(xt02(:,1), Ef2(ntime+1:end)-1.96*sqrt(Varf2(ntime+1:end)), '--', 'color', col2, 'linewidth', 2)
xlabel('age')
title('effect of age for both sexes')


%- Age as baseline hazard

% Age in the beginning and in the end of follow-up
y=[x0(:,1) x0(:,1)+leukemiadata(:,1)/365];
y2=y;
% normalise ages
y2=y2-min(y2(:)); y2=y2./max(y(:));

x2=x;
x2(:,1)=[];

% covariance for proportional part
gpcf = gpcf_sexp('lengthScale', ones(1,size(x2,2)), 'magnSigma2', .5, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% NOTE! if multiple covariance functions per latent is used, define
% gp.comp_cf as follows:
% gp = gp_set(..., 'comp_cf', {[1 2] [5 6]});
% where [1 2] are for hazard function, and [5 6] for proportional part
gp = gp_set('lik', lik, 'cf', {gpcfh gpcf}, 'jitterSigma2', 1e-6, 'comp_cf', {[1] [2]});

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');

opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gp=gp_optim(gp,x2,y2,'z',ye,'opt',opt);

[Ef1, Covf1] = gp_pred(gp, x2, y2, x2, 'z', ye);
Varf1 = diag(Covf1);
xtmpl=linspace(min(y(:)),max(y(:)),50);

figure, hold on, set(gcf, 'color', 'w'),
plot(xtmpl, Ef1(1:ntime),'k', 'linewidth', 3)
plot(xtmpl, Ef1(1:ntime)+1.96*sqrt(Varf1(1:ntime)),'--k', 'linewidth', 2)
plot(xtmpl, Ef1(1:ntime)-1.96*sqrt(Varf1(1:ntime)),'--k', 'linewidth', 2)
title('log-baseline hazard (age)')
xlabel('age')
