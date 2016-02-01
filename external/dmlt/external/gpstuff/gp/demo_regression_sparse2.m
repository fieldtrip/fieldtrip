%DEMO_REGRESSION_SPARSE2  Regression demo comparing different sparse
%                         approximations with optimization of inducing
%                         variables
%
%  Description
%    A regression problem with one input variable and one output
%    variable with Gaussian noise. The output is assumed to be
%    realization of additive functions and Gaussian noise.
% 
%    For standard full GP demonstration, see for example
%    DEMO_REGRESSION1, and for detailed discussion, Rasmussen and
%    Williams (2006). For more sparse demonstrations including use
%    of compact support covariance functions see DEMO_REGRESSION2,
%    DEMO_MODELASSESMENT2, and DEMO_SPARSEAPPROX.
% 
%    In this demo, sparse approximations for the full GP model are
%    compared. We use
%      - FIC, fully independent conditional
%      - DTC, deterministic training conditional
%      - VAR, variational approach
%    For illustration purposes the hyperparameters from the full GP
%    are used for the sparse models and only the inducing variables
%    are optimised.
%     
%    For technical details, see Quinonero-Candela and Rasmussen
%    (2005) for the FIC and DTC models and Titsias (2009) for the
%    VAR model.
% 
%    We use a simple one dimensional data set to present the three
%    methods.
% 
%  See also DEMO_REGRESSION1, DEMO_REGRESSION2, DEMO_REGRESSION_SPARSE1
%
%
%  References:
% 
%    Quinonero-Candela, J. and Rasmussen, C. E. (2005). A Unifying
%    View of Sparse Approximate Gaussian Process Regression. Journal
%    of Machine Learning Research.
% 
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
% 
%    Titsias, M. K. (2009). Variational Model Selection for Sparse
%    Gaussian Process Regression. Technical Report, University of
%    Manchester.

% Copyright (c) 2010 Heikki Peura, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Set randomstream for reproducing same results
prevstream=setrandstream();

% Start by creating 1D data
xx=linspace(1,10,901);

% Choose a subset of data so that the data are less dense in the right end.
% xt are the inputs and yt are the outputs, xstar are the values we want to
% predict.
x1=logspace(0,1,100);
x1=round(x1*100)-99;
x=xx(x1)';
y=2*sin(4*x)+0.2*randn(size(x));
xt=[1:0.01:14]';
[n,nin] = size(x);

fprintf('Full GP\n')
% Initialize full GP with a squared exponential component and set
% priors for their parameters.
pl = prior_t('s2', 1);
pm = prior_logunif();
pn = prior_logunif();

gpcfse = gpcf_sexp('lengthScale',0.5,'magnSigma2',1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
lik = lik_gaussian('sigma2', 0.1, 'sigma2_prior', pn);

gp = gp_set('lik', lik, 'cf', gpcfse, 'jitterSigma2', 1e-6);

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4);
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);

[Eft_full, Varft_full] = gp_pred(gp, x, y, xt);
Varft_full = Varft_full + gp.lik.sigma2;

figure
% Blue crosses are the initial inducing input locations, red ones are
% the optimised ones. Black circles represent the distance to the next
% optimized location, with a dashed trendline.');

subplot(2,2,1);hold on;
plot(xt,Eft_full,'k', 'LineWidth', 2)
plot(xt,Eft_full-2.*sqrt(Varft_full),'--','Color',[0 0.5 0])
plot(xt,Eft_full+2.*sqrt(Varft_full),'--','Color',[0 0.5 0])
plot(x,y,'.', 'MarkerSize',7)
ylim([-3 3])
title('FULL GP')

fprintf('FIC GP\n')
% Run FIC approximation for the same data: choose the inducing
% inputs Xu, then proceed with the inference with the optimized
% parameters from the full GP: here, we optimize only the locations
% of the inducing inputs for the FIC model.
Xu=round(10+90*rand(18,1))/10; % Random placement

% Change type to FIC, add inducing inputs, and optimize only inducing inputs
gp_fic = gp_set(gp, 'type','FIC','X_u',Xu,'infer_params','inducing');

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4);
% Optimize with the scaled conjugate gradient method
gp_fic=gp_optim(gp_fic,x,y,'opt',opt);

[Eft_fic, Varft_fic] = gp_pred(gp_fic, x, y, xt);
Varft_fic = Varft_fic + gp_fic.lik.sigma2;


subplot(2,2,2);hold on;
h1=plot(xt,Eft_fic,'k', 'LineWidth', 2);
h2=plot(xt,Eft_fic-2.*sqrt(Varft_fic),'--','Color',[0 0.5 0]);
plot(xt,Eft_fic+2.*sqrt(Varft_fic),'--','Color',[0 0.5 0])
h3=plot(x,y,'.', 'MarkerSize',7);
h4=plot(Xu, -2.8, 'bx', 'MarkerSize', 5, 'LineWidth', 2);
h5=plot(gp_fic.X_u, -3, 'rx', 'MarkerSize', 5, 'LineWidth', 2);
legend([h1  h2 h3 h4(1) h5(1)],'Ef','95% CI','Data','Initial X_u','Optimized X_u')
% plot diff of sorted X_u and regress line for that
%XuSorted=sort(gp_fic.X_u);
%dXuSorted=diff(XuSorted);
%bb=regress(dXuSorted,[ones(size(dXuSorted)) XuSorted(1:end-1)]);
%plotbb=bb(1)+(min(XuSorted):0.1:max(XuSorted))*bb(2);
%plot(XuSorted(1:end-1),dXuSorted,'ko');
%plot(min(XuSorted):0.1:max(XuSorted),plotbb,'k--')
ylim([-3 3])
title('FIC')


fprintf('VAR GP\n')
% Run the VAR model similarly to the FIC model with the same
% starting inducing inputs. The difference in the optimized results
% is notable. The VAR model places the inducing inputs quite evenly
% (slightly increasing as the data becomes more sparse), with
% predictions closely matching the full GP model. The other two
% sparse approximations yield less reliable results.
gp_var = gp_set(gp,'type','VAR','X_u',Xu,'infer_params','inducing');

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4);
% Optimize with the scaled conjugate gradient method
gp_var=gp_optim(gp_var,x,y,'opt',opt);

[Eft_var, Varft_var] = gp_pred(gp_var, x, y, xt);
Varft_var = Varft_var + gp_var.lik.sigma2;



subplot(2,2,4);hold on
plot(xt,Eft_var,'k', 'LineWidth', 2)
plot(xt,Eft_var-2.*sqrt(Varft_var),'--','Color',[0 0.5 0])
plot(xt,Eft_var+2.*sqrt(Varft_var),'--','Color',[0 0.5 0])
plot(x,y,'.', 'MarkerSize',7)
plot(gp_var.X_u, -3, 'rx', 'MarkerSize', 5, 'LineWidth', 2)
plot(Xu, -2.8, 'bx', 'MarkerSize', 5, 'LineWidth', 2)
% plot diff of sorted X_u and regress line for that
%XuSorted=sort(gp_var.X_u);
%dXuSorted=diff(XuSorted);
%bb=regress(dXuSorted,[ones(size(dXuSorted)) XuSorted(1:end-1)]);
%plotbb=bb(1)+(min(XuSorted):0.1:max(XuSorted))*bb(2);
%plot(XuSorted(1:end-1),dXuSorted,'ko');
%plot(min(XuSorted):0.1:max(XuSorted),plotbb,'k--')
ylim([-3 3])
title('VAR')


fprintf('DTC GP\n')
% Run the DTC model similarly to the FIC model with the same starting
% inducing inputs. The difference in the optimized results is notable.
gp_dtc = gp_set(gp,'type','DTC','X_u',Xu,'infer_params','inducing');

% Set the options for the optimization
opt=optimset('TolFun',1e-4,'TolX',1e-4);
% Optimize with the scaled conjugate gradient method
gp_dtc=gp_optim(gp_dtc,x,y,'opt',opt);

[Eft_dtc, Varft_dtc] = gp_pred(gp_dtc, x, y, xt);
Varft_dtc = Varft_dtc + gp_dtc.lik.sigma2;

subplot(2,2,3);hold on
plot(xt,Eft_dtc,'k', 'LineWidth', 2)
plot(xt,Eft_dtc-2.*sqrt(Varft_dtc),'--','Color',[0 0.5 0])
plot(xt,Eft_dtc+2.*sqrt(Varft_dtc),'--','Color',[0 0.5 0])
plot(x,y,'.', 'MarkerSize',7)
plot(gp_dtc.X_u, -3, 'rx', 'MarkerSize', 5, 'LineWidth', 2)
plot(Xu, -2.8, 'bx', 'MarkerSize', 5, 'LineWidth', 2)
% plot diff of sorted X_u and regress line for that
%XuSorted=sort(gp_dtc.X_u);
%dXuSorted=diff(XuSorted);
%bb=regress(dXuSorted,[ones(size(dXuSorted)) XuSorted(1:end-1)]);
%plotbb=bb(1)+(min(XuSorted):0.1:max(XuSorted))*bb(2);
%plot(XuSorted(1:end-1),dXuSorted,'ko');
%plot(min(XuSorted):0.1:max(XuSorted),plotbb,'k--')
ylim([-3 3])
title('DTC')

% Set back initial random stream
setrandstream([],prevstream);
