%DEMO_BINOMIAL_APC  Demonstration for modeling age-period-cohort data
%                   by a binomial model combined with GP prior.
%
%  Description
%    Demonstration of estimating the unknown the population
%    proportion in binomial model from a sequence of
%    success/failure trials. Data consists of observations Y
%    describing the number of successes in a sequence of N iid
%    (Bernoulli) trials, and of explanatory variables X. The
%    binomial model is
%
%      Y_i ~ Binomial(Y_i | N_i, p_i),
%
%    where the parameter p_i represents the proportion of
%    successes. The total number of trials N_i is fixed in the
%    model. A Gaussian process prior is assumed for latent
%    variables f
%
%      f = N(0, K),
%
%    which are linked to the p_i parameter using the logistic
%    transformation:
%       
%      p_i = logit^-1(f_i) = 1/(1+exp(-f_i)).
%
%    The elements of the covariance matrix K are given as K_ij =
%    k(x_i, x_j | th). The function k(x_i, x_j | th) is covariance
%    function and th its parameters. We place a prior for
%    parameters, p(th). The inference is done with Laplace
%    approximation.
%   
%    In this demonstration X is three dimensional, containing
%    inputs for age group, time period and cohort group (birth
%    year) for each data point. Y represents the number of disease
%    cases among N susceptibles. The data is simulated such that
%    each input have additionally there is a interaction effect
%    between age group and time period.
%
%    NOTE! In the prediction, the total number of trials Nt at the
%    test points Xt must be set additionally in the likelihood
%    structure when E(Yt), Var(Yt) or predictive densities p(Yt)
%    are computed.

% Copyright (c) 2010 Jaakko Riihim�ki, Jouni Hartikainen, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


%========================================================
% data analysis with full GP model
%========================================================

% Demonstration for modeling age-period-cohort data
% by a binomial model combined with GP prior.

% First load data
S = which('demo_binomial_apc');
L = strrep(S,'demo_binomial_apc.m','demodata/binodata.txt');
binodata=load(L);

f = binodata(:,1);
f1 = binodata(:,2);
f2 = binodata(:,3);
f3 = binodata(:,4);
f4 = binodata(:,5);
nn = binodata(:,6);
xx = binodata(:,7:9);
yy = binodata(:,10);

% xx contains the three dimensional inputs for 1377 data points:
%   xx(:,1) - age group
%   xx(:,2) - time period
%   xx(:,3) - cohort (birth year)

% Use only a (random) proportion of original data points in training
inds = randperm(size(xx,1));
ntr = 300;
itr = sort(inds(1:ntr));
itst = sort(inds(ntr+1:end));

% Original points
xxo = xx; yyo = yy;
nno = nn; fo  = f;
f1o = f1; f2o = f2;
f3o = f3; f4o = f4;

% Test points
xt = xx(itst,:); yt = yy(itst,:);
nt = nn(itst,:); ft = f(itst,:);
f1t = f1(itst,:); f2t = f2(itst,:);
f3t = f3(itst,:); f4t = f4(itst,:);

% Training points
xx = xx(itr,:); yy = yy(itr,:);
nn = nn(itr,:); f  = f(itr,:);
f1 = f1(itr,:); f2 = f2(itr,:);
f3 = f3(itr,:); f4 = f4(itr,:);


% Initialization of covariance functions. We shall assume here
% that inputs have additive effect to disease risk (covariance
% functions gpcf1, gpcf2 and gpcf3) as well as interaction effect
% between age group and time period (gpcf4).

% First define priors for length scales and magnitudes
pl = prior_t('s2',10);
pm = prior_sqrtt();
% covariance functions
gpcf1 = gpcf_sexp('selectedVariables', 1,'lengthScale',[10], ...
                  'lengthScale_prior', pl, 'magnSigma2', 1, ...
                  'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('selectedVariables', 2,'lengthScale',[10], ...
                  'lengthScale_prior', pl, 'magnSigma2', 1, ...
                  'magnSigma2_prior', pm);
gpcf3 = gpcf_sexp('selectedVariables', 3,'lengthScale',[4], ...
                  'lengthScale_prior', pl, 'magnSigma2', 1, ...
                  'magnSigma2_prior', pm);
gpcf4 = gpcf_sexp('selectedVariables', [1 2],'lengthScale',[10 2], ...
                  'lengthScale_prior', pl, 'magnSigma2', 1, ...
                  'magnSigma2_prior', pm);

% Initialize the likelihood structure
lik = lik_binomial;
    
% Initialize GP structure
gp = gp_set('lik', lik, 'cf', {gpcf1,gpcf2,gpcf3,gpcf4}, 'jitterSigma2', 1e-6);
    
% Set the approximate inference method
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-2,'TolX',1e-3,'Display','iter');
% Optimize with the BFGS quasi-Newton method
gp=gp_optim(gp,xx,yy,'z',nn,'opt',opt,'optimf',@fminlbfgs);

% Making predictions

% First with all components
[Eft,Varft,lpyt,Eyt,Varyt] = gp_pred(gp,xx,yy,xt,'z',nn,'zt',nt,'yt',yt);
% Age group effect
[Eft_1,Varft_1] = gp_pred(gp,xx,yy,xxo,'predcf',[1],'z',nn,'zt',nno);
% Time period effect
[Eft_2,Varft_2] = gp_pred(gp,xx,yy,xxo,'predcf',[2],'z',nn,'zt',nno);
% Cohort effect
[Eft_3,Varft_3] = gp_pred(gp,xx,yy,xxo,'predcf',[3],'z',nn,'zt',nno);
% Interaction effect between age group and time period
[Eft_4,Varft_4] = gp_pred(gp,xx,yy,xxo,'predcf',[4],'z',nn,'zt',nno);

% Plotting predictions

% First some indexes needed for plotting
% the additive effect
[xx1 ind1] = unique(xxo(:,1));
[xx2 ind2] = unique(xxo(:,2));
[xx3 ind3] = unique(xxo(:,3));

% Age group effect
figure; subplot(3,1,1)
set(gcf, 'color', 'w'), hold on
color1=ones(1,3)*0.8; color2=ones(1,3)*0.5;
% Estimate
h1=fill([xx1' fliplr(xx1')], [(Eft_1(ind1)+1.96*sqrt(Varft_1(ind1)))' ...
    fliplr((Eft_1(ind1)-1.96*sqrt(Varft_1(ind1)))')], color1, 'edgecolor', color1);
h2=plot(xx1, Eft_1(ind1), 'color', color2, 'linewidth', 3);
% True function
h4=plot(xx1, f1o(ind1), 'color', 'r', 'linewidth', 2); hold off
title('Age group effect')
xlabel('Age group'); ylabel('logit(p)')
legend([h4 h2],'True','Estimated')

% Time period effect
subplot(3,1,2)
set(gcf, 'color', 'w'), hold on
h1=fill([xx2' fliplr(xx2')], [(Eft_2(ind2)+1.96*sqrt(Varft_2(ind2)))' ...
    fliplr((Eft_2(ind2)-1.96*sqrt(Varft_2(ind2)))')], color1, 'edgecolor', color1);
h2=plot(xx2, Eft_2(ind2), 'color', color2, 'linewidth', 3);
% true function
h4=plot(xx2, f2o(ind2), 'color', 'r', 'linewidth', 2);
title('Time period effect')
xlabel('Time'); ylabel('logit(p)')
legend([h4 h2],'True','Estimated')

% Cohort effect
subplot(3,1,3)
set(gcf, 'color', 'w'), hold on
h1=fill([xx3' fliplr(xx3')], [(Eft_3(ind3)+1.96*sqrt(Varft_3(ind3)))' ...
    fliplr((Eft_3(ind3)-1.96*sqrt(Varft_3(ind3)))')], color1, 'edgecolor', color1);
h2=plot(xx3, Eft_3(ind3), 'color', color2, 'linewidth', 3);
% true function
h4=plot(xx3, f3o(ind3), 'color', 'r', 'linewidth', 2);
title('Cohort effect')
xlabel('Cohort'); ylabel('logit(p)')
legend([h4 h2],'True','Estimated')

% Plotting of interaction effect
figure; subplot(1,2,1)
imagesc(reshape(Eft_4,81,17)); colorbar;
title('Estimated interaction effect')
xlabel('Time period'); ylabel('Age group')
subplot(1,2,2);
imagesc(reshape(f4o,81,17)); colorbar;
title('True interaction')
xlabel('Time period'); ylabel('Age group')

