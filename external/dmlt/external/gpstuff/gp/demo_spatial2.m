%DEMO_SPATIAL2  Demonstration for a disease mapping problem with Gaussian
%               process prior and negative binomial observation
%               model
%
%  Description
%    The disease mapping problem consist of a data with number of
%    death cases, Y, and background population, N, appointed to
%    co-ordinates, X. The goal is to find a relative risk surface,
%    which describes if the number of death cases in certain areas
%    is lower or higher than expected. The data is simulated.
%
%    The model is constructed as follows:
%
%    The number of death cases Y_i in area i is assumed to satisfy
%
%         Y_i ~ Neg-Bin(Y_i| d, E_i * r_i)
%
%    where E_i is the expected number of deaths (see Vanhatalo and
%    Vehtari (2007), how E_i is evaluated) at area i, r_i is the
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
%    The inference is conducted first with Laplace approximation
%    and then with EP. We use compactly supported covariance
%    function which leads to sparse covariance matrix.
%
%  See also  
%    DEMO_REGRESSION1, DEMO_CLASSIFIC1, DEMO_SPATIAL1

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% =====================================
% Laplace approximation
% =====================================

% load the data
S = which('demo_spatial2');
data = load(strrep(S,'demo_spatial2.m','demodata/spatial2.txt'));

x = data(:,1:2);
ye = data(:,3);
y = data(:,4);

% Now we have loaded the following parameters
% x = co-ordinates 
% y = number of deaths
% ye = the expexted number of deaths

fprintf(['GP with negative-binomial observation model, Laplace\n' ...
         'integration over the latent values and MAP estimate\n' ...
         'for the parameters\n']);

% Create the covariance functions
pl = prior_t();
pm = prior_sqrtt('s2', 0.3);
if ~exist('ldlchol')
  warning('GPstuff:SuiteSparseMissing',...
    ['SuiteSparse is not properly installed. \n' ...
    'Using gpcf_sexp (non-compact support) instead of gpcf_ppcs2 (compact support)']);
  gpcf1 = gpcf_sexp('lengthScale', 5, 'magnSigma2', 0.05);
  gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
else
  gpcf1 = gpcf_ppcs2('nin', 2, 'lengthScale', 5, 'magnSigma2', 0.05);
  gpcf1 = gpcf_ppcs2(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
end

% Create the likelihood structure
lik = lik_negbin();

% Create the GP structure
gp = gp_set('lik', lik, 'cf', gpcf1, 'jitterSigma2', 1e-4); 

% Set the approximate inference method to Laplace
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the scaled conjugate optimization
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt);

% Visualize sparsity pattern
figure
C = gp_trcov(gp,x);
fprintf('Proportion of non-zeros is %.4f\n',nnz(C) / prod(size(C)))
p = amd(C);
spy(C(p,p))

% Make prediction to the data points
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye);

% Define help parameters for plotting
xii=sub2ind([120 70],x(:,2),x(:,1));
[X1,X2]=meshgrid(1:70,1:120);

% Plot the figures
figure
G=repmat(NaN,size(X1));
G(xii)=exp(Ef);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.2    1.5])
axis equal
axis([0 70 0 120])
title('Posterior median of the relative risk (Laplace)')

figure
G=repmat(NaN,size(X1));
G(xii)=(exp(Varf) - 1).*exp(2*Ef+Varf);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
%set(gca, 'Clim', [0.005    0.03])
axis equal
axis([0 70 0 120])
title('Posterior variance of the relative risk (Laplace)')


% =====================================
% EP approximation
% =====================================
fprintf(['GP with negative-binomial observation model, EP\n' ...
         'integration over the latent values and MAP estimate\n' ...
         'for the parameters\n']);

% Set the approximate inference method to EP
gp = gp_set(gp, 'latent_method', 'EP');

% Set the options for the scaled conjugate optimization
opt=optimset('TolFun',1e-2,'TolX',1e-2,'Display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'z',ye,'opt',opt);

% Visualize sparsity pattern
figure
C = gp_trcov(gp,x);
fprintf('Proportion of non-zeros is %.4f\n',nnz(C) / prod(size(C)))
p = amd(C);
spy(C(p,p))

% make prediction to the data points
[Ef, Varf] = gp_pred(gp, x, y, x, 'z', ye);

% Define help parameters for plotting
xii=sub2ind([120 70],x(:,2),x(:,1));
[X1,X2]=meshgrid(1:70,1:120);

% Plot the figures
figure
G=repmat(NaN,size(X1));
G(xii)=exp(Ef);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
set(gca, 'Clim', [0.2    1.5])
axis equal
axis([0 70 0 120])
title('Posterior median of the relative risk (EP)')

figure
G=repmat(NaN,size(X1));
G(xii)=(exp(Varf) - 1).*exp(2*Ef+Varf);
pcolor(X1,X2,G),shading flat
colormap(mapcolor(G)),colorbar
%set(gca, 'Clim', [0.005    0.03])
axis equal
axis([0 70 0 120])
title('Posterior variance of the relative risk (EP)')
