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

% Create the covariance functions
% pl = prior_gaussian('s2',2);
% pm = prior_gaussian('s2',2);
% gpcf1 = gpcf_neuralnetwork('weightSigma2',1, 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);
% gpcf2 = gpcf_neuralnetwork('weightSigma2',1, 'biasSigma2', 0.05, 'weightSigma2_prior', pl, 'biasSigma2_prior', pm);
gpcf1 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, 'lengthScale_prior', prior_logunif(), 'magnSigma2_prior', prior_logunif());
gpcf2 = gpcf_sexp('lengthScale', 1, 'magnSigma2', 1, 'lengthScale_prior', prior_logunif(), 'magnSigma2_prior', prior_logunif());

% Create the likelihood structure
lik = lik_inputdependentweibull('shape', 0.1, 'shape_prior', prior_logunif());

% Create the GP structure
gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'comp_cf', {[1] [2]}, 'jitterSigma2', 1e-6);

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4,'Derivativecheck','off', 'Display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt);

xt1=zeros(200,nin); xt1(:,2)=1;
xt2=zeros(200,nin); xt2(:,2)=-1;

xt1(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt2(:,1)=linspace(min(x(:,1)), max(x(:,1)), 200);
xt01(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);
xt02(:,1)=linspace(min(x0(:,1)), max(x0(:,1)), 200);

% Do predictions
[Ef1, Varf1] = gp_pred(gp, x, y, xt1, 'z', ye);
[Ef2, Varf2] = gp_pred(gp, x, y, xt2, 'z', ye);

% Create normal weibull model
lik = lik_weibull();
gp2 = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);
gp2 = gp_set(gp2, 'latent_method', 'Laplace');
gp2=gp_optim(gp2,x,y,'z',ye,'opt',opt);

[Ef1_2, Varf1_2] = gp_pred(gp2, x, y, xt1, 'z', ye);
[Ef2_2, Varf2_2] = gp_pred(gp2, x, y, xt2, 'z', ye);

lik = lik_loggaussian();
gp3 = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-6);
gp3 = gp_set(gp3, 'latent_method', 'Laplace');
gp3=gp_optim(gp3,x,y,'z',ye,'opt',opt);

[Ef1_3, Varf1_3] = gp_pred(gp3, x, y, xt1, 'z', ye);
[Ef2_3, Varf2_3] = gp_pred(gp3, x, y, xt2, 'z', ye);


% Plot results
col1=ones(1,3)*0.7;
col2=ones(1,3)*0.3;
figure, hold on, set(gcf, 'color', 'w'),
plot(xt01(:,1), Ef1(1:nt), 'color', col1, 'linewidth', 3)
plot(xt01(:,1), Ef1(1:nt)+1.96*sqrt(diag(Varf1(1:nt,1:nt))), '--', 'color', col1, 'linewidth', 2)
plot(xt01(:,1), Ef1(1:nt)-1.96*sqrt(diag(Varf1(1:nt,1:nt))), '--', 'color', col1, 'linewidth', 2)

plot(xt02(:,1), Ef2(1:nt), 'color', col2, 'linewidth', 3)
plot(xt02(:,1), Ef2(1:nt)+1.96*sqrt(diag(Varf2(1:nt,1:nt))), '--', 'color', col2, 'linewidth', 2)
plot(xt02(:,1), Ef2(1:nt)-1.96*sqrt(diag(Varf2(1:nt,1:nt))), '--', 'color', col2, 'linewidth', 2)
xlabel('age')
title('effect of age for both sexes, inputdependent-weibull')

col1=ones(1,3)*0.7;
col2=ones(1,3)*0.3;
figure, hold on, set(gcf, 'color', 'w'),
plot(xt01(:,1), Ef1_2, 'color', col1, 'linewidth', 3)
plot(xt01(:,1), Ef1_2+1.96*sqrt(Varf1_2), '--', 'color', col1, 'linewidth', 2)
plot(xt01(:,1), Ef1_2-1.96*sqrt(Varf1_2), '--', 'color', col1, 'linewidth', 2)

plot(xt02(:,1), Ef2_2, 'color', col2, 'linewidth', 3)
plot(xt02(:,1), Ef2_2+1.96*sqrt(Varf2_2), '--', 'color', col2, 'linewidth', 2)
plot(xt02(:,1), Ef2_2-1.96*sqrt(Varf2_2), '--', 'color', col2, 'linewidth', 2)
xlabel('age')
title('effect of age for both sexes, weibull')

col1=ones(1,3)*0.7;
col2=ones(1,3)*0.3;
figure, hold on, set(gcf, 'color', 'w'),
plot(xt01(:,1), Ef1_3, 'color', col1, 'linewidth', 3)
plot(xt01(:,1), Ef1_3+1.96*sqrt(Varf1_3), '--', 'color', col1, 'linewidth', 2)
plot(xt01(:,1), Ef1_3-1.96*sqrt(Varf1_3), '--', 'color', col1, 'linewidth', 2)

plot(xt02(:,1), Ef2_3, 'color', col2, 'linewidth', 3)
plot(xt02(:,1), Ef2_3+1.96*sqrt(Varf2_3), '--', 'color', col2, 'linewidth', 2)
plot(xt02(:,1), Ef2_3-1.96*sqrt(Varf2_3), '--', 'color', col2, 'linewidth', 2)
xlabel('age')
title('effect of age for both sexes, loggaussian')
