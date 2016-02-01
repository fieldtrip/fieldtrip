%DEMO_BINOMIAL2  Demonstration of Gaussian process model with binomial
%                likelihood
%
%  Description
%    Demonstration of modeling the dose-response relation in
%    bioassay experiment using binomial model with latent linear
%    model. Data consists of observations Y describing the number
%    of successes in a sequence of N iid trials, and log-dose X. 
%    The binomial model is
%
%      Y_i ~ Binomial(Y_i | N_i, p_i),
%
%    where the parameter p_i represents the proportion of
%    successes. The total number of trials N_i is fixed in the
%    model. A Gaussian process prior is assumed for latent
%    variables f (in this example constant+linear model is used)
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
%    NOTE! In the prediction, the total number of trials Nt at the
%    test points Xt must be set additionally in the likelihood
%    structure when E(Yt), Var(Yt) or predictive densities p(Yt)
%    are computed.
%
%    Laplace, EP, and MCMC approximations are compared. Laplace and
%    MCMC results are similar to results in Gelman et al (2004),
%    who used equivalent but explicitly parametrized model.
%
%  Reference
%    Gelman et al (2004). Bayesian data Analysis, second edition,
%      Chapman & Hall/CRC.
%
%  See also DEMO_BINOMIAL1, DEMO_BINOMIAL_APC

% Copyright (c) 2010-2011 Jaakko Riihim�ki, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


% Bioassay data from Gelman et al (p. 89, 2004)
x=[-0.86 -0.3 -0.05 .73]';
% number of trials
N=[5 5 5 5]';
% number of successes
y=[0 1 3 5]';
% n
[n, nin] = size(x);

% equally spaced test points for visualisation
xgrid=linspace(-1.5,1.5,100)';
Ntgrid=ones(size(xgrid))*5;

% Create parts of the covariance function
% Half-Student's t-prior is used to give flat prior in the
% interesting region with respect to the likelihood
cfc = gpcf_constant('constSigma2_prior',prior_t('s2',20^2));
cfl = gpcf_linear('coeffSigma2_prior',prior_t('s2',20^2));
% Create the GP structure
gp = gp_set('lik', lik_binomial(), 'cf', {cfc cfl}, 'jitterSigma2', 1e-8);

% ------- Laplace approximation --------
tic
fprintf('Laplace+grid approximation ')
% (Laplace is default, so this could be skipped)
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'z',N,'opt',opt);
% Form a weighted grid of samples, which is used to integrate over
% the hyperparameters
gpia=gp_ia(gp,x,y,'z',N,'int_method','grid');

% First make predictions just with the linear part to get
% probability that linear effect is positive
[Eft, Varft] = gp_pred(gpia, x, y, xgrid(1), 'z', N, 'zt', Ntgrid(1),'predcf',2);
pbetapositive=normcdf(0,Eft,sqrt(Varft));

% Get samples from the joint distribution of the latent values at 0 and 1
% to compute the corresponding linear model parameters alpha and beta
% in Gelman et al (2004)
fs = gp_rnd(gpia, x, y, [0 1]', 'z', N, 'zt', [5 5]', 'nsamp', 10000);
a=fs(1,:);b=fs(2,:)-fs(1,:);
% compute samples from the LD50 given b>0 (see, Gelman et al (2004))
ld50s=-a(b>0)./b(b>0);

% Make predictions for the latent value at test points
[Eft, Varft] = gp_pred(gpia, x, y, xgrid, 'z', N, 'zt', Ntgrid);

% Visualise the predictions
figure(1)
% a4r
set(gcf,'DefaultAxesPosition',[0.07 0.07 0.88 0.88])
set(gcf, 'color', 'w')
color1=[0.9 0.9 1]; color2=[0 0 0.9];
clf

% Latent function
subplot('Position',[0.03 0.7 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [(Eft+1.96*sqrt(Varft))' fliplr((Eft-1.96*sqrt(Varft))')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, Eft, 'color', color2, 'linewidth', 3);
axis([-1.5 1.5 -14 14])
title('Laplace+grid approximation')
legend([h1 h2],'Mean of the latent','95% CI',2)
legend('boxoff')
line(xlim,[0 0],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
xlabel('log dose')

% Probability and observations
subplot('Position',[0.03 0.38 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [(logitinv(Eft+1.96*sqrt(Varft)))' fliplr((logitinv(Eft-1.96*sqrt(Varft)))')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, logitinv(Eft), 'color', color2, 'linewidth', 3);
% observations
h3=plot(x, y./5, 'xr', 'markersize', 10, 'linewidth', 2);
legend([h1 h2 h3],'Expected prob.','95% CI','Observations',2)
legend('boxoff')
line(xlim,[.5 .5],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
axis([-1.5 1.5 0 1])
xlabel('log dose')

% Histogram of LD50
subplot('Position',[0.03 0.06 0.3 0.25])
hist(ld50s,-0.6:.04:0.6),set(gca,'xlim',[-.6 .6])
h=get(gca,'Children');
% set(h,'FaceColor',color1)
set(gca,'ytick',[])
ylim([0 2000])
xlabel('LD50')
h1=text(prctile(ld50s,2.5),1850,'2.5%','HorizontalAlignment','center');
h2=text(prctile(ld50s,50),1850,'50%','HorizontalAlignment','center');
h3=text(prctile(ld50s,97.5),1850,'97.5%','HorizontalAlignment','center');
hl=line(repmat(prctile(ld50s,[2.5 50 97.5]),2,1),repmat([0 1800]',1,3),'Color','k');
fprintf('Elapsed time %.0fs\n',toc)

% ------- EP approximation --------
tic
fprintf('EP+grid approximation      ')

gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'z',N,'opt',opt);
% Form a weighted grid of samples, which is used to integrate over
% the hyperparameters
gpia=gp_ia(gp,x,y,'z',N,'int_method','grid');

% First make predictions just with the linear part to get
% probability that linear effect is positive
[Eft, Varft] = gp_pred(gpia, x, y, xgrid(1), 'z', N, 'zt', Ntgrid(1),'predcf',2);
pbetapositive=normcdf(0,Eft,sqrt(Varft));

% Get samples from the joint distribution of the latent values at 0 and 1
% to compute the corresponding linear model parameters alpha and beta
% in Gelman et al (2004)
fs = gp_rnd(gpia, x, y, [0 1]', 'z', N, 'zt', [5 5]', 'nsamp', 10000);
a=fs(1,:);b=fs(2,:)-fs(1,:);
% compute samples from the LD50 given b>0 (see, Gelman et al (2004))
ld50s=-a(b>0)./b(b>0);

% Make predictions for the latent value at test points
[Eft, Varft] = gp_pred(gpia, x, y, xgrid, 'z', N, 'zt', Ntgrid);

% Visualise the predictions

% Latent function
subplot('Position',[0.36 0.7 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [(Eft+1.96*sqrt(Varft))' fliplr((Eft-1.96*sqrt(Varft))')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, Eft, 'color', color2, 'linewidth', 3);
axis([-1.5 1.5 -14 14])
title('EP+grid approximation')
legend([h1 h2],'Mean of the latent','95% CI',2)
legend('boxoff')
line(xlim,[0 0],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
xlabel('log dose')

% Probability and observations
subplot('Position',[0.36 0.38 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [(logitinv(Eft+1.96*sqrt(Varft)))' fliplr((logitinv(Eft-1.96*sqrt(Varft)))')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, logitinv(Eft), 'color', color2, 'linewidth', 3);
% observations
h3=plot(x, y./5, 'xr', 'markersize', 10, 'linewidth', 2);
legend([h1 h2 h3],'Expected prob.','95% CI','Observations',2)
legend('boxoff')
line(xlim,[.5 .5],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
axis([-1.5 1.5 0 1])
xlabel('log dose')

% Histogram of LD50
subplot('Position',[0.36 0.06 0.3 0.25])
hist(ld50s,-0.6:.04:0.6),set(gca,'xlim',[-.6 .6])
h=get(gca,'Children');
% set(h,'FaceColor',color1)
set(gca,'ytick',[])
ylim([0 2000])
xlabel('LD50')
h1=text(prctile(ld50s,2.5),1850,'2.5%','HorizontalAlignment','center');
h2=text(prctile(ld50s,50),1850,'50%','HorizontalAlignment','center');
h3=text(prctile(ld50s,97.5),1850,'97.5%','HorizontalAlignment','center');
hl=line(repmat(prctile(ld50s,[2.5 50 97.5]),2,1),repmat([0 1800]',1,3),'Color','k');
fprintf('Elapsed time %.0fs\n',toc)

% ------- MCMC approximation --------
tic
fprintf('MCMC approximation         ')

gp = gp_set(gp, 'latent_method', 'MCMC', 'jitterSigma2', 1e-4);

[rgp,g,opt] = gp_mc(gp, x, y, 'z', N, 'nsamples', 500, 'repeat', 4, 'display', 0);
rgp=thin(rgp,101);

% First make predictions just with the linear part to get
% probability that linear effect is positive
[Efts, Varfts] = gpmc_preds(rgp, x, y, xgrid(1), 'z', N, 'zt', Ntgrid(1),'predcf',2);
pbetapositive=mean((Efts+randn(size(Efts)).*sqrt(Varfts))<0);

% Get samples from the joint distribution of the latent values at 0 and 1
% to compute the corresponding linear model parameters alpha and beta
% in Gelman et al (2004)
fs = gp_rnd(rgp, x, y, [0 1]', 'z', N, 'zt', [5 5]', 'nsamp', 10000);
a=fs(1,:);b=fs(2,:)-fs(1,:);
% compute samples from the LD50 given b>0 (see, Gelman et al (2004))
ld50s=-a(b>0)./b(b>0);

% Make predictions for the latent value at test points
[Efts, Varfts] = gpmc_preds(rgp, x, y, xgrid, 'z', N, 'zt', Ntgrid);
fs=Efts+randn(size(Efts)).*sqrt(Varfts);

% Visualise the predictions

% Latent function
subplot('Position',[0.69 0.7 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [prctile(fs,97.5,2)' fliplr(prctile(fs,2.5,2)')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, mean(fs,2), 'color', color2, 'linewidth', 3);
axis([-1.5 1.5 -14 14])
title('MCMC approximation')
legend([h1 h2],'Mean of the latent','95% CI',2)
legend('boxoff')
line(xlim,[0 0],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
xlabel('log dose')

% Probability and observations
subplot('Position',[0.69 0.38 0.3 0.25])
hold on
% GP 95% credible interval
h2=fill([xgrid' fliplr(xgrid')], [(logitinv(prctile(fs,97.5,2)))' fliplr(logitinv(prctile(fs,2.5,2))')], color1, 'edgecolor', color1);
% GP mean
h1=plot(xgrid, logitinv(mean(fs,2)), 'color', color2, 'linewidth', 3);
% observations
h3=plot(x, y./5, 'xr', 'markersize', 10, 'linewidth', 2);
legend([h1 h2 h3],'Expected prob.','95% CI','Observations',2)
legend('boxoff')
line(xlim,[.5 .5],'Linestyle','--','color','k')
line([0 0],ylim,'Linestyle','--','color','k')
axis([-1.5 1.5 0 1])
xlabel('log dose')

% Histogram of LD50
subplot('Position',[0.69 0.06 0.3 0.25])
hist(ld50s,-0.6:.04:0.6),set(gca,'xlim',[-.6 .6]);
h=get(gca,'Children');
% set(h,'FaceColor',color1)
set(gca,'ytick',[])
ylim([0 2000])
xlabel('LD50')
h1=text(prctile(ld50s,2.5),1850,'2.5%','HorizontalAlignment','center');
h2=text(prctile(ld50s,50),1850,'50%','HorizontalAlignment','center');
h3=text(prctile(ld50s,97.5),1850,'97.5%','HorizontalAlignment','center');
hl=line(repmat(prctile(ld50s,[2.5 50 97.5]),2,1),repmat([0 1800]',1,3),'Color','k');
fprintf('Elapsed time %.0fs\n',toc)
