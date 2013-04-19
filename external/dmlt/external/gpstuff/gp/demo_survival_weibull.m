%DEMO_SURVIVAL_WEIBULL  Survival model using Weibull baseline hazard
%
%  Description 
%    Survival model using Weibull distribution for the baseline hazard. The
%    hazard rate is 
%   
%       h(t) = h_0(t)*exp(f),
%
%    where the baseline hazard is assumed to be a Weibull distribution
%
%       h_0(t) = r*t^(r-1), r>0
%
%    where r is the shape parameter. A zero-mean Gaussian process prior is
%    assumed for for f = [f_1, f_2,...,f_n] ~ N(0, K), where K is the
%    covariance matrix, whose elements are given as K_ij = k(x_i, x_j |
%    th). The function k(x_i, x_j | th) is covariance function and th its
%    parameters.
%
%    The inference is conducted via EP or Laplace, where we find
%    Gaussian approximation for p(f| th, data), where th is the
%    maximum a posterior (MAP) estimate for the parameters.
%    
%    Example data set is leukemia survival data in Northwest England
%    presented in (Henderson, R., Shimakura, S., and Gorst, D. (2002).
%    Modeling spatial variation in leukemia survival data. Journal of the
%    American Statistical Association, 97:965–972). Data set was downloaded
%    from http://www.math.ntnu.no/%7Ehrue/r-inla.org/examples/leukemia/leuk.dat
%
%  See also  DEMO_SPATIAL1

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
                        %                  ye = 1 for right censored event

% choose (for example) 'age', 'sex', 'wbc', and 'tpi' covariates
x0=leukemiadata(:,5:8);
x=x0;
% normalize continuous covariates 
x(:,[1 3:4])=bsxfun(@rdivide,bsxfun(@minus,x0(:,[1 3:4]),mean(x0(:,[1 3:4]),1)),std(x0(:,[1 3:4]),1));

[n, nin]=size(x);

% Create the covariance functions
pl = prior_gaussian('s2',1);
pm = prior_gaussian('s2',1);
gpcf1 = gpcf_neuralnetwork('weightSigma2', 1*ones(1,nin), 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);

% Create the likelihood structure
lik = lik_weibull();
%lik = lik_loglogistic();
%lik = lik_loggaussian();

% Create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');
%gp = gp_set(gp, 'latent_method', 'EP');

%gradcheck(gp_pak(gp),@gpla_e,@gpla_g,gp,x,y,'z',ye)
%gradcheck(gp_pak(gp),@gpep_e,@gpep_g,gp,x,y,'z',ye)

% Set the options for the optimization
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt,'optimf',@fminlbfgs);

% Make prediction
xt1=zeros(200,nin); xt1(:,2)=1;
xt2=zeros(200,nin); xt2(:,2)=-1;

xt1(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt2(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt01(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);
xt02(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);

[Ef1, Varf1] = gp_pred(gp, x, y, xt1, 'z', ye);
[Ef2, Varf2] = gp_pred(gp, x, y, xt2, 'z', ye);

% Compute posterior mean and 95% credible intervals of latent function as a
% function of age and sex
col1=ones(1,3)*0.7;
col2=ones(1,3)*0.3;
figure, hold on, set(gcf, 'color', 'w'),
plot(xt01(:,1), exp(-Ef1), 'color', col1, 'linewidth', 3)
plot(xt01(:,1), exp(-Ef1+1.96*sqrt(Varf1)), '--', 'color', col1, 'linewidth', 2)
plot(xt01(:,1), exp(-Ef1-1.96*sqrt(Varf1)), '--', 'color', col1, 'linewidth', 2)

plot(xt02(:,1), exp(-Ef2), 'color', col2, 'linewidth', 3)
plot(xt02(:,1), exp(-Ef2+1.96*sqrt(Varf2)), '--', 'color', col2, 'linewidth', 2)
plot(xt02(:,1), exp(-Ef2-1.96*sqrt(Varf2)), '--', 'color', col2, 'linewidth', 2)
set(gca,'YTick',[1:15])
xlabel('age')
title('effect of age for both sexes')

