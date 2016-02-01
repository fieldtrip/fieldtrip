%  DEMO_MODEL_COMPARISON  Model Assessment and Comparisons  
%
%  Description: 
%    
%    we compare the predictive ability of different models by
%    estimating various assessment statistics such as Harrel's C or
%    IDI. Predictions are made using k-fold-CV and uncertainty
%    using Bayesian bootstrap.
%
%    We will compare four different models: Cox proportional hazards,
%    Weibull, log-Gaussian and log-logistic model.   
%   
%    The censoring indicator ye is
%    
%      ye = 0 for uncensored event
%      ye = 1 for right censored event.
% 
%    Example data set is leukemia survival data in Northwest
%    England presented in (Henderson, R., Shimakura, S., and Gorst,
%    D. (2002). Modeling spatial variation in leukemia survival
%    data. Journal of the American Statistical Association,
%    97:965-972). Data set was downloaded from
%    http://www.math.ntnu.no/%7Ehrue/r-inla.org/examples/leukemia/leuk.dat
%
%  See also  DEMO_MODELCOMPARISON, DEMO_SURVIVAL_COXPH
%

% Copyright (c) 2012 Ernesto Ulloa
% Copyright (c) 2012 Aki Vehtari

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

%  we choose for the new model: 'age', 'sex', 'wbc', and 'tpi'
%  covariates
x0=leukemiadata(:,5:8);
x=x0;

% normalize continuous covariates
x(:,[1 3:4])=bsxfun(@rdivide,bsxfun(@minus,x0(:,[1 3:4]),mean(x0(:,[1 3:4]),1)),std(x0(:,[1 3:4]),1));
[n, nin]=size(x);

%* set tau
tt=0.1:.1:1;

% set D event indicator vector for each time in tt (Di=0 if i experienced
% the event before tau and Di=1 otherwise

% Also we set YY, the observed time vector for each time value in tt 
for i=1:size(tt,2)
  for i2=1:size(ye,1)
    if y(i2)>tt(i)
      yytemp(i2)=tt(i);
      Dtemp(i2)=1;   
    else
      if ye(i2)==1
        Dtemp(i2)=1;
      else  
        Dtemp(i2)=0;
      end
      yytemp(i2)=y(i2);
    end
  end
  yyi{i}=yytemp';
  Di{i}=Dtemp';
end
for i=1:size(Di,2)
  D(:,i)=Di{i};
end
for i=1:size(yyi,2)
  yy(:,i)=yyi{i};
end

% set time vector to make predictions
yt=bsxfun(@times,ones(size(y)),tt);

% *** Cox proportional Harzards model ***

% number of time intervals for Cox proportional hazards model
ntime=50;
% create finite partition of time axis
S=linspace(0,max(y)+0.001,ntime+1);

% Create the covariance functions
pl = prior_t('s2',1, 'nu', 4);
pm = prior_t('s2',1, 'nu', 4); 

% covariance for hazard function
gpcfh1 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcfh2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% covariance for proportional part
gpcf1 = gpcf_sexp('lengthScale', ones(1,size(x,2)), 'magnSigma2', 1.2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('lengthScale', ones(1,size(x,2)), 'magnSigma2', 1.2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
lik = lik_coxph('S', S);

% NOTE! if Multiple covariance functions per latent is used, define
% gp.comp_cf as follows:
% gp.comp_cf = {[1 2] [5 6]}
% where [1 2] are for hazard function, and [5 6] for proportional part
gpcph = gp_set('lik', lik, 'cf', {gpcfh1 gpcf2}, 'jitterSigma2', 1e-6);
gpcph.comp_cf = {[1] [2]};

% Set the approximate inference method to Laplace
gpcph = gp_set(gpcph, 'latent_method', 'Laplace');

opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter','Derivativecheck','off');

%  First obtain predictions using k-fold-CV
[cdf]=gp_kfcv_cdf(gpcph,x,y,'z',D,'yt',yt,'opt',opt);
critcph=cdf;

% *** Weibull model ***

% Create the likelihood structure
lik = lik_weibull();

% Create the covariance functions
pl = prior_gaussian('s2',1);
pm = prior_gaussian('s2',1);
gpcf1 = gpcf_neuralnetwork('weightSigma2', 1*ones(1,nin), 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);

% Create the GP structure
gpw = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);

% Set the approximate inference method to Laplace
gpw = gp_set(gpw, 'latent_method', 'Laplace');

% Obtain predictions
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter','Derivativecheck','off');

[cdf]=gp_kfcv_cdf(gpw,x,y,'z',D,'yt',yt,'opt',opt);
critw=cdf;

% *** Log Gaussian model ***

% Create the likelihood structure
lik = lik_loggaussian();

% Create the covariance functions
pl = prior_gaussian('s2',1);
pm = prior_gaussian('s2',1);
gpcf1 = gpcf_neuralnetwork('weightSigma2', 1*ones(1,nin), 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);

% Create the GP structure
gplg = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);

% Set the approximate inference method to Laplace
gplg = gp_set(gplg, 'latent_method', 'Laplace');

% Obtain predictions
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter','Derivativecheck','off');

[cdf]=gp_kfcv_cdf(gplg,x,y,'z',D,'yt',yt,'opt',opt);
critlg=cdf;

% *** Log Logistic model ***

% Create the likelihood structure
lik = lik_loglogistic();

% Create the covariance functions
pl = prior_gaussian('s2',1);
pm = prior_gaussian('s2',1);
gpcf1 = gpcf_neuralnetwork('weightSigma2', 1*ones(1,nin), 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);

% Create the GP structure
gpll = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);

% Set the approximate inference method to Laplace
gpll = gp_set(gpll, 'latent_method', 'Laplace');

% Obtain predictions
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter','Derivativecheck','off');

[cdf]=gp_kfcv_cdf(gpll,x,y,'z',D,'yt',yt,'opt',opt);
critll=cdf;

%% Calculate statics and compare models 
% MODEL COMPARISON

%% AUC 
% AUC for Binary outcomes P(Pi>Pj | Di=1,Dj=0)
[auccph,fpscph,tpscph]=aucs(critcph(:,length(tt)),D(:,length(tt)));
[aucw,fpsw,tpsw]=aucs(critw(:,length(tt)),D(:,length(tt)));
[aucll,fpsll,tpsll]=aucs(critll(:,length(tt)),D(:,length(tt)));
[auclg,fpslg,tpslg]=aucs(critlg(:,length(tt)),D(:,length(tt)));

fprintf(['\n AUC at end of study for Cox-ph:   ', num2str(auccph)]);
fprintf(['\n AUC at end of study for Weibull:   ', num2str(aucw)]);
fprintf(['\n AUC at end of study for log-Gaussian:   ', num2str(auclg)]);
fprintf(['\n AUC at end of study for log-logistic:   ', num2str(aucll)]);

figure
hold on
subplot(2,2,1);
plot(fpscph,tpscph,'b')
legend('Cox Ph')
subplot(2,2,2);
plot(fpsw,tpsw,'r')
legend('Weibull')
subplot(2,2,3);
plot(fpsll,tpsll,'g')
legend('Log Logistic')
subplot(2,2,4);
plot(fpslg,tpslg,'y')
title('ROC curve')
legend('Log Gaussian')
hold off

%% Harrell's C 
%Harrel's C(t) = P(Pi>Pj | Di(ti)=1, ti<tj, ti<tt) at end of study
ccph=hct(critcph(:,length(tt)),y,D(:,length(tt)),tt(length(tt)));
cw=hct(critw(:,length(tt)),y,D(:,length(tt)),tt(length(tt)));
cll=hct(critll(:,length(tt)),y,D(:,length(tt)),tt(length(tt)));
clg=hct(critlg(:,length(tt)),y,D(:,length(tt)),tt(length(tt)));

fprintf(['\n Harrells C at end of study for Cox-ph:   ', num2str(ccph)]);
fprintf(['\n Harrells C at end of study for Weibull:   ', num2str(cw)]);
fprintf(['\n Harrells C at end of study for log-Gaussian:   ', num2str(cll)]);
fprintf(['\n Harrells C at end of study for log-logistic:   ', num2str(clg)]);

% Obtain for Log Gaussian and Log Logistic models Binary AUC(t) = P(Pi>Pj | Di(t)=1,Dj(t)=0) and
[autlgll,clgll]=assess(critlg,critll,yy,D,tt);

% Again, but for Cox-ph and Weibull models
[autcphw,ccphw]=assess(critcph,critw,yy,D,tt);

% Plot Harrells C in function of time for all four models, time interval:(0,1) 
plot(tt,clgll(:,1),'r');
hold on;
plot(tt,clgll(:,2),'g');
plot(tt,ccphw(:,1),'b');
plot(tt,ccphw(:,2),'y');
legend('Log Gaussian','Log Logistic','Cox Proportional Hazards Model','Weibull')
title('Harrells C in function of time ');
xlabel('Time');
ylabel('Harrolds C');
hold off;

% Use Bayesian bootsrap to obtain Harrells (CLG-CLL) statistic
% density at tt=1
[clg,bbll]=hcs(critll(:,size(tt,2)),y,D(:,size(tt,2)),tt(size(tt,2)),'rsubstream',1);
[cll,bblg]=hcs(critlg(:,size(tt,2)),y,D(:,size(tt,2)),tt(size(tt,2)),'rsubstream',1);
title('Estimated density of CLG-CLL')
hold on 
lgpdens(bblg-bbll)
hold off 

% We integrate the (CLG-CLL) estimated density in the (0,inf) interval
zc=lgpdens_cum(bblg-bbll,0,inf);
fprintf(['Estimated c statistics for Log Gaussian and Log Logistic respectively:   ', num2str(clg) '  ' num2str(cll)]);
fprintf(['Cumulative probability in the (0,inf) interval:   ', num2str(zc)]);

% Again, use Bayesian bootsrap to obtain Harrells (CW-CCPH)
% statistic density at tt=1
[cw,bbw]=hcs(critw(:,size(tt,2)),y,D(:,size(tt,2)),tt(size(tt,2)),'rsubstream',1);
[ccph,bbcph]=hcs(critcph(:,size(tt,2)),y,D(:,size(tt,2)),tt(size(tt,2)),'rsubstream',1);
title('Estimated density of CCW-CCPH')
hold on 
lgpdens(bbw-bbcph)
hold off 

% We integrate the (CW-CCPH) estimated density in the (0,inf) interval
zc=lgpdens_cum(bbw-bbcph,0,inf);
fprintf(['Estimated C-statistics for Weibull and Cox-ph respectively:   ', num2str(cw) '  ' num2str(ccph)]);
fprintf(['Cumulative probability in the (0,inf) interval:   ', num2str(zc)]);

%% IDI

% Estimate R^2 for all four models, We calculate IDI between log
% logistic and loggaussian and IDI between weibull and coxph. Also
% we estimate the densities and the cumulative probability in the
% (0,inf) interval at time 1

[idi1,bbid1,rll,rlg] = idis(critll(:,size(tt,2)),critlg(:,size(tt,2)),'rsubstream',1);
zidi1=lgpdens_cum(bbid1,0,inf);

fprintf(['\n R^2 statistic for log-logistic model:', num2str(rll)]);
fprintf(['\n R^2 statistic for log-Gaussian model:', num2str(rlg)]);
display(stll)
display(stlg)

fprintf(['Estimated IDI between log-logistic and log-Gaussian : ', num2str(idi1)]);
fprintf(['cumulative probability in the (0,inf) interval: ', num2str(zidi1)]);
display(st1)
display(st2)

title('IDI estimated density between log-logistic and log-Gaussian')
hold on 
lgpdens(bbid1)
hold off 

[idi2,bbid2,rw,rcph] = idis(critw(:,size(tt,2)),critcph(:,size(tt,2)),'rsubstream',1);
zidi2=lgpdens_cum(bbid2,0,inf);


fprintf(['\n R^2 statistic for Weibull model:', num2str(rw)]);
fprintf(['\n R^2 statistic for Cox-ph model:', num2str(rcph)]);
display(st1)
display(st2)

fprintf(['Estimated idi between Weibull and Cox-ph : ', num2str(idi2)]);
fprintf(['cumulative probability in the (0,inf) interval: ', num2str(zidi2)]);
display(st1)
display(st2)

title('IDI estimated density between Weibull and Cox-ph')
hold on 
lgpdens(bbid2)
hold off 

%% Ext AUC

% Ext_AUC for different subsets of tt 
Indxtmp{1}=1:1:size(tt,2);
Indx{1}=1:1:size(tt,2);
j=2;
k=round(size(tt,2)/2); 
for i=2:k
  Indxtmp{i}=1:i:size(tt,2);
  if length(Indxtmp{i})~=length(Indxtmp{i-1})
    Indx{j}=Indxtmp{i};
    j=j+1;
  end
  
end

for i=1:size(Indx,2)
  l(i)=length(Indx{i});
end

for i=1:size(Indx,2)
  eacph(i) = ext_auc(critcph(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
  eaw(i) = ext_auc(critw(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
  ealg(i) = ext_auc(critlg(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
  eall(i) = ext_auc(critll(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
end

hold on
xlabel('Number of distinct time partitions')
ylabel('Extended AUC')
plot(wrev(l),eacph,'r')
plot(wrev(l),eaw,'b')
plot(wrev(l),ealg,'g')
plot(wrev(l),eall,'y')
legend('Cox-ph', 'Weibull','Log-Gaussian','Log-Logistic')
hold off

stcph=sprintf(['\n ExtAUC at end of study for model 1:   ', num2str(eacph(size(Indx,2)))]);
stw=sprintf(['\n ExtAUC at end of study for model 2:   ', num2str(eaw(size(Indx,2)))]);
stll=sprintf(['\n ExtAUC at end of study for model 1:   ', num2str(eall(size(Indx,2)))]);
stlg=sprintf(['\n ExtAUC at end of study for model 2:   ', num2str(ealg(size(Indx,2)))]);

display(stcph)
display(stw)
display(stll)
display(stlg)

%% plot predictions for average individual 
%********************  Superimpose a plot of a prediction for an average
% individual
% choose (for example) 'age', 'sex', 'wbc', and 'tpi' covariates
% average covariates except sex 
xa=mean(leukemiadata(:,[5,7,8]),1);
% -1 for female
xa=[xa(1) -1 xa(2:3)];
% [Ef1, Varf1] = gp_pred(gp, x, y,xa ,'z', ye);
% 1 for male 
xa=[xa(1) 1 xa(2:3)];

% *** Obtain predictions 


% Cox Proportional Hazards model

% optimise parameters
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gpcph=gp_optim(gpcph,x,y,'z',ye,'opt',opt);
[hcph survcph]=pred_coxphhs(gpcph,x,y,xa,'z',ye);
yycph=[1 survcph];


% Weibull model

% optimise parameters
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gpw=gp_optim(gpw,x,y,'z',ye,'opt',opt);

% obtain samples for average individual 
fs_w = gp_rnd(gpw, x, y, xa , 'z', ye,'nsamp',1000);
zz=[gpcph.lik.stime];
%zz = 0.001:.001:max(y);
yyw = mean(exp(-bsxfun(@times,exp(fs_w'),zz.^(gpw.lik.shape))));


% Log Gaussian model 

% optimize parameters 
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gplg=gp_optim(gplg,x,y,'z',ye,'opt',opt);

% obtain samples for average individual 
fs_lg = gp_rnd(gplg, x, y, xa , 'z', ye,'nsamp',1000);

%zz = 0.001:.001:max(y);
for i=1:size(zz,2)   
  yylg(i) =1. - mean(norm_cdf(log(zz(i)), fs_lg, sqrt(gplg.lik.sigma2)),2);
end

% Log Logistic model 

% optimize parameters 
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
gpll=gp_optim(gpll,x,y,'z',ye,'opt',opt);

% obtain samples for average individual 
fs_ll = gp_rnd(gpll, x, y, xa , 'z', ye,'nsamp',1000);

%zz = 0.001:.001:max(y);
for i=1:size(zz,2)
  yyll(i)=1.-mean(1./(1.+exp(-(-(fs_ll)+log(zz(i)))/gpll.lik.shape)),2);
end
%yyll=(1+(zz/exp(mean(fs_ll))).^gpll.lik.shape).^(-1);




% Calculate and plot survival function of empirical and the four different models 
[f,z] = ecdf(y,'censoring',ye,'function','survivor');
clf
stairs(z,f,'k','LineWidth',2)
hold on
% stairs(z,flo,'r:','LineWidth',2)
% stairs(z,fup,'r:','LineWidth',2)
zzcph=[gpcph.lik.stime];
plot(zz,yycph,'r--','LineWidth',2)
plot(zz,yyw,'g--','LineWidth',2)
plot(zz,yylg,'b--','LineWidth',2)
plot(zz,yyll,'y--','LineWidth',2)
legend('Empirical','Coxph','Weibull','Log-gaussian','Log-logistic')
title('Survival predictions for average individual ');
hold off





