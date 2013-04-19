%DEMO_SURVIVAL_COMPARISON  Survival model comparison
%
%  Description: 
%    
%    By using kfc-validation and Bayesian bootstrap we compare the
%    predictive ability of different models by estimating various
%    assessment statistics .
%
%    We will compare two Cox proportional hazars model, the first
%    model will have less covariates than the second model.
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
%    97:965–972). Data set was downloaded from
%    http://www.math.ntnu.no/%7Ehrue/r-inla.org/examples/leukemia/leuk.dat
%
%  See also  DEMO_SURVIVAL_COMPARISON2, DEMO_SURVIVAL_COXPH

% Copyright (C) 2012 Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

%% First load data
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

%  we choose for the first model: 'age' and 'sex'covariates
x01=leukemiadata(:,5:6);
x1=x01;

%  we choose for the second model: 'age', 'sex', 'wbc', and 'tpi' covariates
x02=leukemiadata(:,5:8);
x2=x02;

% normalize continuous covariates 

x1(:,1)=normdata(x01(:,1));
x2(:,[1 3:4])=normdata(x02(:,[1 3:4]));

[n1, nin1]=size(x1);
[n2, nin2]=size(x2);

% number of time intervals
ntime=50;
% create finite partition of time axis
S=linspace(0,max(y)+0.001,ntime+1);

%% obtain predictions

% Create the covariance functions
pl = prior_t('s2',1, 'nu', 4);
pm = prior_t('s2',1, 'nu', 4); 

% covariance for hazard function
gpcfh1 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcfh2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% covariance for proportional part
gpcf1 = gpcf_sexp('lengthScale', ones(1,size(x1,2)), 'magnSigma2', 1.2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gpcf2 = gpcf_sexp('lengthScale', ones(1,size(x2,2)), 'magnSigma2', 1.2, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
lik = lik_coxph('S', S);

gp1 = gp_set('lik', lik, 'cf', {gpcfh1 gpcf1}, 'jitterSigma2', 1e-6, 'comp_cf', {[1] [2]});
gp2 = gp_set('lik', lik, 'cf', {gpcfh2 gpcf2}, 'jitterSigma2', 1e-6, 'comp_cf', {[1] [2]});

% Set the approximate inference method to Laplace
gp1 = gp_set(gp1, 'latent_method', 'Laplace');
gp2 = gp_set(gp2, 'latent_method', 'Laplace');

opt=optimset('TolFun',1e-2,'TolX',1e-4,'Display','iter','Derivativecheck','off');

% obtain predictions for both models using kfc-validation

%* first we set tau
tt=0.1:.1:1;

% set D event indicator vector for each time in tt (Di=0 if i experienced
% the event before tau and Di=1 otherwise)
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

% Obtain predictions
% (This takes several minutes)
crit1=gp_kfcv_cdf(gp1,x1,y,'z',D,'yt',yt,'opt',opt);
crit2=gp_kfcv_cdf(gp2,x2,y,'z',D,'yt',yt,'opt',opt);

%% Calculate statics and compare models 
% MODEL COMPARISON

%% AUC 
% AUC for Binary outcomes P(Pi>Pj | Di=1,Dj=0)
[auc1,fps1,tps1]=aucs(crit1(:,length(tt)),D(:,length(tt)));
[auc2,fps2,tps2]=aucs(crit2(:,length(tt)),D(:,length(tt)));

fprintf('AUC at end of study for model 1: %.3f \n', auc1);
fprintf('AUC at end of study for model 2: %.3f \n', auc2);
hold on
plot(fps1,tps1,'b')
plot(fps2,tps2,'r')
title('ROC curve')
legend('model 1', 'model 2',4)
xlabel('False positives')
ylabel('True positives')
hold off

%% Harrell's C
% Obtain for both models Binary AUC(t) = P(Pi>Pj | Di(t)=1,Dj(t)=0) and
% Harrell's C(t) = P(Pi>Pj | Di(ti)=1, ti<tj, ti<tt) for every element of tt  
ct1=hct(crit1,yy,D,tt);
ct2=hct(crit2,yy,D,tt);
auct1=auct(crit1,yy,D,tt);
auct2=auct(crit2,yy,D,tt);
c=[ct1 ct2];

% Plot for both models Harrells C in function of time
plot(tt,c(:,1),'r');
hold on;
plot(tt,c(:,2),'g');
legend('Old model','New model')
title('Harrolds C in function of time ');
xlabel('Time');
ylabel('Harrell''s C');
hold off;

%% Estimated density
% Use bayesian bootsrap to obtain Harrells (C1-C2) statistic density at tt=1
[c1,bb1]=hcs(crit1(:,end),y,ye,1,'rsubstream',1);
[c2,bb2]=hcs(crit2(:,end),y,ye,1,'rsubstream',1);
title('Estimated density of C2-C1')
lgpdens(bb2-bb1)
xlabel('Difference in Harrell''s C statistics (C2-C1)');

% We integrate the (C1-C2) estimated density in the (0,inf) interval
zc=lgpdens_cum(bb2-bb1,0,inf);
fprintf('Estimated c statistics for model 1 and 2 respectively:  %.3f, %.3f \n', c1, c2);
fprintf('cumulative probability in the (0,inf) interval:  %.2f \n', zc);

%% IDI
%Estimate R^2 for both models, idi, its density and the cumulative
%probability in the (0,inf) interval, al at time 1

[idi,r1,r2,bbid] = idis(crit1(:,end),crit2(:,end),'rsubstream',1);
zidi=lgpdens_cum(bbid,0,inf);
title('IDI estimated density')
lgpdens(bbid)

fprintf('R^2 statistic for model 1: %.3f \n', r1);
fprintf('R^2 statistic for model 2: %.3f \n', r2);

fprintf('Estimated idi: %.3f ', idi);
fprintf('cumulative probability in the (0,inf) interval: %.2f\n', zidi);

%% EXT AUC

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
  ea1(i) = ext_auc(crit1(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
  ea2(i) = ext_auc(crit2(:,Indx{i}),tt(:,Indx{i}),tt(:,Indx{i}(size(Indx{i},2))));
end

extauc1 = ext_auc(crit1,tt,tt(:,end));
extauc2 = ext_auc(crit2,tt,tt(:,end));

fprintf('ExtAUC at end of study for model 1: %.3f \n', extauc1);
fprintf('ExtAUC at end of study for model 2: %.3f \n', extauc2);
