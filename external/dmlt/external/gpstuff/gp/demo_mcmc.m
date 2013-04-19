%DEMO_MCMC  Demonstration of Markov Chain Monte Carlo sampling methods in
%           GPstuff
%
%  Description
%    This demo consists of two parts. First part demonstrates sampling of
%    hyperparameters and using Laplace approximation for latent inference,
%    while the second part demonstrates sampling of hyperparameters and
%    latent values. The data is drowning deaths in Finland  and we use
%    Poisson likelihood (see demo_periodic for more thorough explanation).
%
%    Sampling methods in this demo include HMC, HMC-NUTS, Slice
%    sampling, Adaptive slice sampling (for the first part) and surrogate 
%    slice sampling (for the second part). All other samplers besides
%    surrogate slice sampler work in both just sampling hyperparameters and
%    sampling hyperparameters with latent values.
%
%    Note that this demo is not intended for comparing 
%    how efficient different samplers are (different option parameters are 
%    not optimized for this problem), but rather as an example as in how 
%    to use them.
%
%  See also 
%    DEMO_PERIODIC, GP_MC, HMC2, HMC_NUTS, SLS, SURROGATE_SLS
%

% Copyright (c) 2012 Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

S = which('demo_periodic');
L = strrep(S,'demo_periodic.m','demodata/drowning.txt');
data=load(L);
y = data(:, 2:13);
y=y';
y=y(:);
y1=y;
y=y(1:72);
avgy = mean(y);
x = [1:length(y)]';

[n,nin] = size(x);

% --- Construct the model ---

gpcf1 = gpcf_sexp('lengthScale', [67], 'magnSigma2', 1);
gpcfp = gpcf_periodic('lengthScale', [1.3], 'magnSigma2', 2.4*2.4,...
    'period', 12,'lengthScale_sexp', 50, 'decay', 1);
likn=gpcf_neuralnetwork('biasSigma2',10,'weightSigma2',3);
gpcf2 = gpcf_sexp('lengthScale', [2], 'magnSigma2', 2);

% ... Then set the prior for the parameters of covariance functions...
pl = prior_t('s2', 1000, 'nu', 3);
pm = prior_sqrtt('s2', 2, 'nu', 3);
pl2 = prior_t('s2', 5, 'nu', 3);
pm2 = prior_sqrtt('s2', 3, 'nu', 3);
ppl = prior_t('s2', 100, 'nu', 3);
ppm = prior_sqrtt('s2', 1, 'nu', 3);
pn = prior_t('s2', 10, 'nu', 4);
ppp = prior_t('s2', 100, 'nu', 4);

gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% gpcf2 = gpcf_sexp(gpcf2, 'lengthScale_prior', pl2, 'magnSigma2_prior', pm2);
likn = gpcf_neuralnetwork(likn, 'biasSigma2_prior', pn, 'weightSigma2_prior', ppp);
gpcfp = gpcf_periodic(gpcfp, 'lengthScale_prior', ppl, 'magnSigma2_prior', ppm,  'lengthScale_sexp_prior', pl);

% ... Create the GP structure, Poisson likelihood
z=repmat(mean(y),length(y),1);
gp = gp_set('lik', lik_poisson(), 'cf', {gpcf1,gpcfp,likn}, 'jitterSigma2', 1e-6);

opt=optimset('TolX',1e-4,'TolFun',1e-4); 
gp=gp_optim(gp,x,y,'z',z,'opt',opt);

%========================================================
% PART 1 Sampling of hyperparameters
%========================================================

% Here we set latent method to Laplace (this could be omitted as Laplace is
% default latent method)
gp = gp_set(gp, 'latent_method', 'Laplace');

% ---------------------------------
% Hybrid Monte Carlo (HMC) sampling
% ---------------------------------
fprintf('Sampling of hyperparameters with HMC.\n')


% Set some options for HMC
hmc_opt.nsamples=1;
hmc_opt.decay=0.8;
hmc_opt.persistence=0;
hmc_opt.stepadj=0.04;
hmc_opt.steps=10;

% gp_mc does the actual sampling. 'Repeat' (default 1) is the number of samples we
% generate before accepting one, 'hmc_opt', defines the option structure
% and also sets the sampler for HMC/NUTS and 'nsamples' is the number of samples
% we would like for the function to return. Note that real number of
% samples generated depends on sampler, 'repeat' and 'nsamples'.
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'repeat',1, 'hmc_opt', hmc_opt, 'nsamples', 100);

% Here we remove the burn-in in the samples
rgp=thin(rgp, 20);

% Here we save the acquired hyperparameter samples
w1_hmc=gp_pak(rgp);

% clear options structure
clear('hmc_opt');

% ----------------------------
% No-U-Turn Sampler (HMC-NUTS)
% ----------------------------
fprintf('Sampling of hyperparameters with HMC-NUTS.\n')

% For NUTS, we only need to specify that we use NUTS-HMC instead of
% original HMC, and the number of adaptation steps of step-size parameter
hmc_opt.nuts = 1; % True or false
hmc_opt.nadapt = 20; % Number of step-size adaptation steps (Burn-in)

% Generate samples
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'repeat',1, 'hmc_opt', hmc_opt, 'nsamples', 100);

% Options structure opt includes the adapted step-size parameters opt.hmc_opt.epsilon
% and more diagnostics for evaluating convergence of the samples.

rgp=thin(rgp, 20);

% Save samples
w1_nuts=gp_pak(rgp);

clear('opt','hmc_opt');

% --------------
% Slice sampling
% --------------
fprintf('Sampling of hyperparameters with Slice sampler.\n')

% Set options for slice sampler (for more thorough explanation, see
% sls_opt.m)
sls_opt.nomit = 0;
sls_opt.display = 0;

% Use 'stepping' method
% sls_opt.method = 'stepping';
opt.method = 'minmax';
% opt.method = 'multimm';

sls_opt.overrelaxation = 0;
sls_opt.alimit = 4;
sls_opt.wsize = 1;
sls_opt.mlimit = 4;
sls_opt.maxiter = 50;
sls_opt.plimit = 2;
sls_opt.unimodal = 0;
sls_opt.mmlimits = [sls_opt.wsize-(sls_opt.wsize*sls_opt.mlimit); sls_opt.wsize+(sls_opt.wsize*sls_opt.mlimit)];

% Do the actual sampling
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'repeat',1, 'sls_opt', sls_opt, 'nsamples', 100);
% Note that since slice sampler is default sampler for sampling ONLY
% hyperparameters, we could use easier command below
% [rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'repeat',2, 'nsamples', 100);

rgp=thin(rgp, 20);

% Save hyperparameters
w1_sls=gp_pak(rgp);

clear('opt','sls_opt')

% -----------------------
% Adaptice slice sampling
% -----------------------
fprintf('Sampling of hyperparameters with Shrinking Rank Slice sampler.\n')

% Set options for adaptive slice samplers (see sls_opt.m)
opt.method='shrnk';
% opt.method='covmatch';
opt.sigma = 1;
opt.display = 0;

% Do the actual sampling
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'repeat',1, 'sls_opt', sls_opt, 'nsamples', 100);

rgp=thin(rgp, 20);

% Save hyperparameters
w1_asls=gp_pak(rgp);

clear('opt','sls_opt', 'latent_opt');

%========================================================
% PART 2 Sampling of hyperparameters and latent values
%========================================================

% Here we set latent method to MCMC
gp = gp_set(gp, 'latent_method', 'MCMC');

% Sampler for latent values can be set with
latent_opt.method=@esls;
% latent_opt.method=@scaled_hmc;
% latent_opt.method=@scaled_mh;
gp = gp_set(gp, 'latent_opt', latent_opt);
% Here we use elliptical slice sampler (esls), which is also default
% sampler for latent values (meaning that we could omit setting the
% sampler)
clear('latent_opt');

% ---------------------------------
% Hybrid Monte Carlo (HMC) sampling
% ---------------------------------
fprintf('Sampling of hyperparameters with HMC and latent values with Elliptical slice sampler.\n')
% Options for HMC
hmc_opt.nsamples=1;
hmc_opt.decay=0.8;
hmc_opt.persistence=0;
hmc_opt.stepadj=0.03;
hmc_opt.steps=20;

% As we sample also latent values now, we must set options for that sampler
% as well if we do not want to use defaulsavet options

% esls only needs the repeat option, which we set to 30 here. This is 
% equivalent for running with repeat = 1 and discarding every sample 
% besides every 10th
latent_opt.repeat = 30;

% Do the sampling, latent_opt should be given if user doesn't want to use
% default options for default latent sampler.
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'hmc_opt', hmc_opt, 'latent_opt', latent_opt, 'nsamples', 100);

rgp=thin(rgp, 20);

% Here we save the acquired hyperparameter samples
w2_hmc=gp_pak(rgp);

% clear options structure
clear('hmc_opt');

% ----------------------------
% No-U-Turn Sampler (HMC-NUTS)
% ----------------------------
fprintf('Sampling of hyperparameters with HMC-NUTS and latent values with Elliptical slice sampler.\n')

% Options for NUTS
hmc_opt.nuts = 1; 
hmc_opt.nadapt = 20; 

% Use same latent options

% Generate samples
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'hmc_opt', hmc_opt, 'latent_opt', latent_opt, 'nsamples', 100);

rgp=thin(rgp, 20);

% Save samples
w2_nuts=gp_pak(rgp);

clear('opt','hmc_opt', 'latent_opt');

% --------------
% Slice sampling
% --------------
fprintf('Sampling of hyperparameters with Slice sampler and latent values with HMC.\n')

% Options for slice sampler
sls_opt.nomit = 0;
sls_opt.display = 0;

% sls_opt.method = 'stepping';
sls_opt.method = 'minmax';
% sls_opt.method = 'multimm';

sls_opt.overrelaxation = 0;
sls_opt.alimit = 4;
sls_opt.wsize = 1;
sls_opt.mlimit = 4;
sls_opt.maxiter = 50;
sls_opt.plimit = 2;
sls_opt.unimodal = 0;
sls_opt.mmlimits = [sls_opt.wsize-(sls_opt.wsize*sls_opt.mlimit); sls_opt.wsize+(sls_opt.wsize*sls_opt.mlimit)];

% Here we demonstrate different latent value sampler, scaled_hmc
gp=gp_set(gp,'latent_opt', struct('method', @scaled_hmc));

% HMC-latent options
latent_opt.nsamples=1;
latent_opt.nomit=0;
latent_opt.persistence=0;
latent_opt.window=5;
latent_opt.repeat=1;
latent_opt.steps=7;
latent_opt.window=1;
latent_opt.stepadj=0.15;

% Do the actual sampling
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'sls_opt', sls_opt, 'latent_opt', latent_opt, 'nsamples', 100);

rgp=thin(rgp, 20);

% Save hyperparameters
w2_sls=gp_pak(rgp);

clear('opt','sls_opt')

% ------------------------
% Surrogate slice sampling
% ------------------------
fprintf('Sampling of hyperparameters with Surrogate slice sampler and latent values with Elliptical slice sampler.\n')

% Use elliptical slice sampler
gp = gp_set(gp, 'latent_opt', struct('method', @esls));

% Options for surrogate slice sampler (same options as in sls)
ssls_opt.method = 'stepping';
% sls_opt.method = 'minmax';
% sls_opt.method = 'multimm';

ssls_opt.overrelaxation = 0;
ssls_opt.alimit = 4;
ssls_opt.wsize = 1;
ssls_opt.mlimit = 4;
ssls_opt.maxiter = 50;
ssls_opt.plimit = 2;
ssls_opt.unimodal = 0;
ssls_opt.mmlimits = [ssls_opt.wsize-(ssls_opt.wsize*ssls_opt.mlimit); ssls_opt.wsize+(ssls_opt.wsize*ssls_opt.mlimit)];

% Surrogate slice sampler does both the sampling of hyperparameters and
% latent values, so we need to give the sampler also options for latent
% value sampling
ssls_opt.latent_opt.repeat=20;

% Do the actual sampling, now we give options as a ssls_opt. Note that we
% dont give latent_opt, because surrogate sls takes care of the latent
% value sampling
[rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'ssls_opt', ssls_opt, 'nsamples', 100);

rgp=thin(rgp, 20);

% As surrogate slice sampler is default sampler when sampling both
% hyperparameters and latent values, we could use the easier command below 
% [rgp,g,opt]=gp_mc(gp, x, y, 'z', z, 'nsamples', 100);

% Save hyperparameters
w2_ssls=gp_pak(rgp);


% Plot some sampling results
figure,
subplot(4,1,1),plot(w1_hmc),  title('Sampling only hyperparameters, HMC'); ylim([-10 10]), xlim([0 80]);
subplot(4,1,2),plot(w1_nuts), title('HMC-NUTS'); ylim([-10 10]), xlim([0 80]);
subplot(4,1,3),plot(w1_sls), title('Slice sampling'); ylim([-10 10]), xlim([0 80]);
subplot(4,1,4), plot(w1_asls), title('Adaptive Slice sampling (Shrinking-Rank)'); ylim([-10 10]), xlim([0 80]);

figure,
subplot(4,1,1),plot(w2_hmc),  title('Sampling hyperparameters and latent values, HMC');ylim([-10 10]), xlim([0 80]);
subplot(4,1,2),plot(w2_nuts), title('HMC-NUTS');ylim([-10 10]), xlim([0 80]);
subplot(4,1,3),plot(w2_sls), title('Slice sampling'); ylim([-10 10]), xlim([0 80]);
subplot(4,1,4),plot(w2_ssls), title('Surrogate Slice sampling'); ylim([-10 10]), xlim([0 80]);

