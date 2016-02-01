% DEMO_LGPDENS  Demonstration of Logistic-Gaussian Process density estimate
%               for 1D and 2D data 
%
%  Description 
%
%    Logistic Gaussian Process (LGPDENS) is a model for density
%    estimation. For the samples from continuous distribution, the
%    space is discretized into n intervals with equal lengths
%    covering the interesting region. The following model is used
%    in estimation
%    
%        p(y_i|f_i) ~ exp(f_i) / Sum_j^n exp(f_j),
%
%    where a zero mean Gaussian process prior is placed for f =
%    [f_1, f_2,...,f_n] ~ N(0, K). K is the covariance matrix,
%    whose elements are given as K_ij = k(x_i, x_j | th). The
%    function k(x_i, x_j| th) is covariance function and th its
%    parameters, hyperparameters. We place a hyperprior for
%    hyperparameters, p(th).
%
%    The inference is conducted via Laplace and the last example
%    compares the results of Laplace approximation to MCMC.
%

% Copyright (c) 2011 Jaakko Riihimäki and Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

prevstream=setrandstream(0,'mrg32k3a');

% =====================================
% 1) 1D-examples
% =====================================

figure(1)
subplot(2,2,1)
% t_4
setrandstream(0,'mrg32k3a');
x=trnd(4,1,100)';
xt=linspace(-7,7,400)';
lgpdens(x,xt);
axis tight
title('t_4')
% true density
p0=t_pdf(xt,4,0,1);
line(xt,p0,'color','k')
%sum(p0.*log(p))

subplot(2,2,2)
% Mixture of two t_4
setrandstream(0,'mrg32k3a');
n1=sum(rand(100,1)<3/4);
n2=100-n1;
x=[trnd(4,n1,1); 3+trnd(4,n2,1)/4];
xt=linspace(-6,6,400)';
lgpdens(x,xt);
axis tight
title('Mixture of two t_4')
% true density
p0=t_pdf(xt,4,0,1)*2/3+t_pdf(xt,4,3,1/4)*1/3;
line(xt,p0,'color','k')

subplot(2,2,3)
% Galaxy data
S = which('demo_lgpdens');
L = strrep(S,'demo_lgpdens.m','demodata/galaxy.txt');
x=load(L);
xt=linspace(0,40000,200)';
lgpdens(x,xt);
axis tight
title('Galaxy data')
% true density is unknown

subplot(2,2,4)
% Gamma(1,1)
setrandstream(0,'mrg32k3a');
x=gamrnd(1,1,100,1);
xt=linspace(0,5,400)';
lgpdens(x,xt);
axis tight
title('Gamma(1,1)')
p0=gam_pdf(xt,1,1);
line(xt,p0,'color','k')


% =====================================
% 1) 2D-examples
% =====================================

figure(2)
clf
subplot(2,2,1)
% t_4
n=100;
Sigma = [1 .7; .7 1];R = chol(Sigma);
setrandstream(0,'mrg32k3a');
x=trnd(8,n,2)*R;
lgpdens(x);
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
axis([-4 4 -4 4])
title('Student t_4')

subplot(2,2,2)
% Old faithful
L = strrep(S,'demo_lgpdens.m','demodata/faithful.txt');
x=load(L);
lgpdens(x,'range',[1 6 40 100]);
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
title('Old faithful')

subplot(2,2,3)
% Banana-shaped
n=100;
setrandstream(0,'mrg32k3a');
b=0.02;x=randn(n,2);x(:,1)=x(:,1)*10;x(:,2)=x(:,2)+b*x(:,1).^2-10*b;
lgpdens(x,'range',[-30 30 -5 20],'gridn',26);
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
axis([-25 25 -5 10])
title('Banana')

subplot(2,2,4)
% Ring
n=100;
setrandstream(0,'mrg32k3a');
phi=(rand(n,1)-0.5)*2*pi;
x=[1.5*cos(phi)+randn(n,1)*0.2 1.5*sin(phi)+randn(n,1)*0.2];
lgpdens(x,'gridn',30);
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
axis([-2.5 2.5 -2.5 2.5])
title('Ring')


% =====================================
% 1) 2D-example of conditional density estimate
% =====================================
figure(3)
% Ring
lgpdens(x,'gridn',30, 'cond_dens', 'on');
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
axis([-2.5 2.5 -2.5 2.5])
title('Ring - conditional density')


% =====================================
% 1D-example with FFT speed-up and 2D-example
% with Kronecker product (low-rank) speed-up
% =====================================
figure(4)
% Galaxy data with finer grid, 1D
L = strrep(S,'demo_lgpdens.m','demodata/galaxy.txt');
x=load(L);
xt=linspace(0,40000,800)';
subplot(2,2,1)
tic,lgpdens(x,xt,'speedup', 'off');t0=toc;
axis tight
title(['Galaxy, no speed-up, elapsed time: ' num2str(t0)])
subplot(2,2,2)
tic,lgpdens(x,xt,'speedup', 'on');t1=toc;
axis tight
title(['Galaxy, FFT speed-up, elapsed time: ' num2str(t1)])

% Old faithful, 2D
L = strrep(S,'demo_lgpdens.m','demodata/faithful.txt');
x=load(L);
subplot(2,2,3)
tic,lgpdens(x,'range',[1 6 40 100],'gridn', 30, 'speedup', 'off');t0=toc;
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
title(['Old faithful, no speed-up, elapsed time: ' num2str(t0)])
subplot(2,2,4)
tic,lgpdens(x,'range',[1 6 40 100],'gridn', 30, 'speedup', 'on');t1=toc;
line(x(:,1),x(:,2),'LineStyle','none','Marker','.')
title(['Old faithful, KRON speed-up, elapsed time: ' num2str(t1)])

% =====================================
% 1) 1D-example MCMC vs Laplace
% =====================================
figure(5)
clf
subplot(2,1,1)
% t_4
setrandstream(0,'mrg32k3a');
x=[trnd(4,1,100)]';
xt=linspace(-6,6,200)';
[p,pq]=lgpdens(x,xt);
pla=p;
line(xt,p,'color','r','marker','none','linewidth',2)
line(xt,pq,'color','r','marker','none','linewidth',1,'linestyle','--')
xlim([-7 7])
title('t_4 (Laplace)')
% true density
p0=t_pdf(xt,4,0,1);
line(xt,p0,'color','k')

subplot(2,1,2)
[p,pq]=lgpdens(x,xt,'latent_method','MCMC');
pmc=p;
line(xt,p,'color','r','marker','none','linewidth',2)
line(xt,pq,'color','r','marker','none','linewidth',1,'linestyle','--')
xlim([-7 7])
title('t_4 (MCMC)')
line(xt,p0,'color','k')

[pks] = ksdensity(x,xt);

disp(['Laplace: ' num2str(sum(p0.*log(pla)))])
disp(['MCMC: ' num2str(sum(p0.*log(pmc)))])
disp(['ksdensity: ' num2str(sum(p0.*log(pks)))])

% Set back initial random stream
setrandstream(prevstream);