%DEMO_SPATIAL1  Demonstration for a disease mapping problem
%               with Gaussian process prior and negative-Binomial
%               likelihood
%
%  Description
%    The disease mapping problem consist of a data with number of
%    death cases, Y, and background population, N, appointed to
%    co-ordinates, X. The goal is to find a relative risk surface,
%    which describes if the number of death cases in certain areas
%    is lower or higher than expected. The data consists of the
%    heart attacks in Finland from 1996-2000 aggregated into
%    20kmx20km lattice cells.
%
%    The model constructed is as follows:
%
%    The number of death cases Y_i in area i is assumed to satisfy
%
%         Y_i ~ Neg-Bin(Y_i| d, E_i * r_i)
%
%    where E_i is the expected number of deaths (see Vanhatalo and
%    Vehtari (2007, 2010), how E_i is evaluated) at area i, r_i is the
%    relative risk and d is the dispersion parameter coverning the
%    variance.
%
%    We place a zero mean Gaussian process prior for log(R), R =
%    [r_1, r_2,...,r_n], which implies that at the observed input
%    locations latent values, f, have prior
%
%         f = log(R) ~ N(0, K),
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters. We place a prior
%    for parameters, p(th).
%
%    Since the data set used in this demo is rather large we use
%    FIC sparse approximation for the GP prior.
%
%    The inference is conducted first with Laplace approximation
%    and then via MCMC. We sample from the full posterior p(f, th|
%    data) by alternating the sampling from the conditional
%    posteriors p(f | th, data) and p(th | f, data). The sampling
%    from the conditional posteriors is done by hybrid Monte Carlo
%    (see, for example, Neal, 1996).
%
%    See Vanhatalo and Vehtari (2007) and Vanhatalo et.al. (2010)
%    for more detailed discussion.
%
%    This demo is organised in three parts:
%     1) data analysis with Laplace approximation
%     2) data analysis with integrated Laplace approximation
%     3) data analysis with MCMC
%
%    See also  DEMO_REGRESSION1, DEMO_CLASSIFIC1
%
%  References:
%    Vanhatalo, J., Pietil�inen V. and Vehtari, A. (2010). 
%    Approximate inference for disease mapping with sparse Gaussian
%    processes. Statistics in Medicine, 29(15):.
%
%    Jarno Vanhatalo and Aki Vehtari (2007). Sparse Log Gaussian
%    Processes via MCMC for Spatial Epidemiology. JMLR Workshop and
%    Conference Proceedings, 1:73-89. (Gaussian Processes in
%    Practice)

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


% =====================================
% 1) Laplace approximation
% =====================================

% load the data
S = which('demo_spatial1');
data = load(strrep(S,'demo_spatial1.m','demodata/spatial1.txt'));

x = data(:,1:2);
ye = data(:,3);
y = data(:,4);

% Now we have loaded the following parameters
% x = co-ordinates 
% y = number of deaths
% ye = the expexted number of deaths

% Set the inducing inputs in a regular grid.  Set_PIC returns the
% inducing inputs and block indeces for PIC. It also plots the data
% points, inducing inputs and blocks.
dims = [1    60     1    35];
[trindex, Xu] = set_PIC(x, dims, 5, 'corners', 0);

[n,nin] = size(x);

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_sqrtunif();
gpcf1 = gpcf_matern32('lengthScale', 1, 'magnSigma2', 0.1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
%gpcf2 = gpcf_ppcs3('nin',nin,'lengthScale', 5, 'magnSigma2', 0.05, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);

% Create the likelihood structure
% The data is overdispersed compared to Poisson-model, and with
% with Poisson model the posterior of the hyperparameters is
% multimodal. Negative-Binomial is more robust alternative.
% lik = lik_poisson(); 
lik = lik_negbin();

% Create the FIC GP structure so that inducing inputs are not optimized
gp = gp_set('type', 'FIC', 'lik', lik, 'cf', gpcf1, 'X_u', Xu, ...
            'jitterSigma2', 1e-4, 'infer_params', 'covariance');
% alternative models
%gp = gp_set('type', 'PIC', 'lik', lik, 'cf', gpcf1, 'X_u', Xu, ...
%            'jitterSigma2', 1e-4, 'infer_params', 'covariance',...
%            'tr_index', trindex);
%gp = gp_set('type', 'CS+FIC', 'lik', lik, 'cf', {gpcf1 gpcf2}, 'X_u', Xu, ...
%            'jitterSigma2', 1e-4, 'infer_params', 'covariance');

% --- MAP estimate with Laplace approximation ---

% Set the approximate inference method to Laplace approximation
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the quasi-Newton optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3);
% Optimize with the quasi-Newton method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt);

% make prediction to the data points
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye, 'tstind', [1:n]);

% Define help parameters for plotting
xii=sub2ind([60 35],x(:,2),x(:,1));
[X1,X2]=meshgrid(1:35,1:60);

% Plot the figures
% In the results it should be noticed that:
% - there is much more people living in the south than in the north. 
%   This results in rather high variance in the north
% - The eastern Finland is known to be worse than western Finland in 
%   heart diseases also from other studies.
% - The inducing inputs are set slightly too sparsely for this data, 
%   which results in oversmoothness in the maps
figure
G=repmat(NaN,size(X1));
G(xii)=exp(Ef);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.8    1.25])
axis equal
axis([0 35 0 60])
title('Posterior median of the relative risk, FIC')

figure
G=repmat(NaN,size(X1));
G(xii)=sqrt((exp(Varf) - 1).*exp(2*Ef+Varf));
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.035    0.125])
axis equal
axis([0 35 0 60])
title('Posterior std of the relative risk, FIC')

% the MAP estimate of the parameters. 
S1=sprintf('MAP: length-scale: %.1f, magnSigma: %.3f \n', gp.cf{1}.lengthScale, sqrt(gp.cf{1}.magnSigma2))

% --- IA ---
[gp_array,pth,th,Ef,Varf] = gp_ia(gp, x, y, x, 'z', ye, ...
                               'int_method', 'grid', 'step_size',2);

% Plot the figures
figure
G=repmat(NaN,size(X1));
G(xii)=exp(Ef);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.8    1.25])
axis equal
axis([0 35 0 60])
title('Posterior median of the relative risk, FIC')

figure
G=repmat(NaN,size(X1));
G(xii)=sqrt((exp(Varf) - 1).*exp(2*Ef+Varf));
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.035    0.125])
axis equal
axis([0 35 0 60])
title('Posterior std of the relative risk, FIC')

% the IA estimate of the parameters 
% in real life
Elth=sum(bsxfun(@times,pth,th));
Elth2=sum(bsxfun(@times,pth,th.^2));
Stdlth = sqrt(Elth2-Elth.^2);
S2=sprintf('IA-GRID: 90%% CI - length-scale: [%.1f,%.1f], magnSigma: [%.3f,%.3f] \n',exp(Elth(2)-1.645*Stdlth(2)),exp(Elth(2)+1.645*Stdlth(2)), sqrt(exp(Elth(1)-1.645*Stdlth(1))),sqrt(exp(Elth(1)+1.645*Stdlth(1))))

% --- MCMC ---

% Set the approximate inference method to MCMC
gp = gp_set(gp, 'latent_method', 'MCMC', 'latent_opt', struct('method',@scaled_hmc));

% Set the sampling options

% HMC-param
hmc_opt.steps=3;
hmc_opt.stepadj=0.01;
hmc_opt.nsamples=1;
hmc_opt.persistence=0;
hmc_opt.decay=0.8;
    
% HMC-latent
latent_opt.nsamples=1;
latent_opt.nomit=0;
latent_opt.persistence=0;
latent_opt.repeat=20;
latent_opt.steps=20;
latent_opt.stepadj=0.15;
latent_opt.window=5;

% Here we make an initialization with 
% slow sampling parameters
[rgp,gp,opt]=gp_mc(gp, x, y, 'z', ye, 'repeat',5, 'hmc_opt', hmc_opt, 'latent_opt', latent_opt);

% Now we reset the sampling parameters to 
% achieve faster sampling
opt.latent_opt.repeat=1;
opt.latent_opt.steps=7;
opt.latent_opt.window=1;
opt.latent_opt.stepadj=0.15;
opt.hmc_opt.persistence=0;
opt.hmc_opt.stepadj=0.03;
opt.hmc_opt.steps=2;

opt.display = 1;
opt.hmc_opt.display = 0;
opt.latent_opt.display=0;

% Define help parameters for plotting
xii=sub2ind([60 35],x(:,2),x(:,1));
[X1,X2]=meshgrid(1:35,1:60);

% Conduct the actual sampling.
% Inside the loop we sample one sample from the latent values and
% parameters at each iteration. After that we plot the samples so
% that we can visually inspect the progress of sampling
figure
% first prepare figures for faster plotting
clf
subplot(1,2,1)
h1=plot(rgp.cf{1}.lengthScale, sqrt(rgp.cf{1}.magnSigma2),...
        rgp.cf{1}.lengthScale(end), sqrt(rgp.cf{1}.magnSigma2(end)),'r*')
xlabel('length-scale')
ylabel('magnitude')
subplot(1,2,2)
G=repmat(NaN,size(X1));
G(xii)=exp(gp.latentValues);
h2=pcolor(X1,X2,G);shading flat
colormap(mapcolor(G)),colorbar
axis equal
axis([0 35 0 60])
title('relative risk')
while length(rgp.edata)<500 %
  [rgp,gp,opt]=gp_mc(gp, x, y, 'record', rgp, 'z', ye, opt, 'display', 0);
  fprintf('        mean hmcrej: %.2f latrej: %.2f\n', mean(rgp.hmcrejects), mean(rgp.lrejects))
  set(h1(1),'XData',rgp.cf{1}.lengthScale,'YData',sqrt(rgp.cf{1}.magnSigma2));
  set(h1(2),'XData',rgp.cf{1}.lengthScale(end),'YData',sqrt(rgp.cf{1}.magnSigma2(end)));
  G=repmat(NaN,size(X1));
  G(xii)=exp(gp.latentValues);
  set(h2,'XData',X1,'YData',X2,'CData',G);
  set(gca,'clim',[0.8 1.25])
  drawnow
end

figure
clf
G=repmat(NaN,size(X1));
G(xii)=median(exp(rgp.latentValues));
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.8    1.25])
axis equal
axis([0 35 0 60])
title('Posterior median of relative risk, FIC GP')

figure
G=repmat(NaN,size(X1));
G(xii)=std(exp(rgp.latentValues), [], 1);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.035    0.125])
axis equal
axis([0 35 0 60])
title('Posterior std of relative risk, FIC GP')

S3=sprintf('MCMC: 90%% CI - length-scale: [%.1f,%.1f], magnSigma: [%.3f,%.3f] \n',prctile(rgp.cf{1}.lengthScale(56:end),5),prctile(rgp.cf{1}.lengthScale(56:end),95),prctile(sqrt(rgp.cf{1}.magnSigma2(56:end)),5),prctile(sqrt(rgp.cf{1}.magnSigma2(56:end)),95))

disp(S1)
disp(S2)
disp(S3)

% $$$ figure
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'DefaultAxesPosition',[0.0  0.02   0.85   0.96]);
% $$$ set(gcf,'DefaultAxesFontSize',8)   %6 8
% $$$ set(gcf,'DefaultTextFontSize',8)   %6 8
% $$$ G=repmat(NaN,size(X1));
% $$$ G(xii)=exp(Ef);
% $$$ pcolor(X1,X2,G),shading flat
% $$$ colormap(mapcolor(G, [1 1])),colorbar
% $$$ axis equal
% $$$ axis([0 35 0 60])
% $$$ set(gca,'YTick',[])
% $$$ set(gca,'XTick',[])
% $$$ set(gca,'XTicklabel',[])
% $$$ set(gca,'YTickLabel',[])
% $$$ set(gcf,'pos',[13.6    10   5.4  6])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
% $$$ 
% $$$ figure
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'DefaultAxesPosition',[0.0  0.02   0.85   0.96]);
% $$$ set(gcf,'DefaultAxesFontSize',8)   %6 8
% $$$ set(gcf,'DefaultTextFontSize',8)   %6 8
% $$$ G=repmat(NaN,size(X1));
% $$$ G(xii)=(exp(Varf) - 1).*exp(2*Ef+Varf);
% $$$ pcolor(X1,X2,G),shading flat
% $$$ colormap(mapcolor(G)),colorbar
% $$$ axis equal
% $$$ axis([0 35 0 60])
% $$$ set(gca,'YTick',[])
% $$$ set(gca,'XTick',[])
% $$$ set(gca,'XTicklabel',[])
% $$$ set(gca,'YTickLabel',[])
% $$$ set(gcf,'pos',[13.6    10   5.5  6])
% $$$ set(gcf,'paperunits',get(gcf,'units'))
% $$$ set(gcf,'paperpos',get(gcf,'pos'))
