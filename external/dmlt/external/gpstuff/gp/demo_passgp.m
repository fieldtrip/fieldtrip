%DEMO_PASSGP Demonstration of PASS-GP routine for GP classification
%
%  Description
%    Here we demonstrate pass-gp routine for Gaussian Processes
%    classification. Data used is 2-dimensional toy data with Gaussian
%    bumbs defining classes. We demonstrate with both fixed and not fixed
%    sizes of active set for pass-gp.
%
%  See also PASSGP

% Copyright (c) 2013 Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Generate toy data
prevstream=setrandstream(0);
[x1,x2]=meshgrid(-5:0.1:5,-5:0.1:5);
x=[x1(:) x2(:)]; x=x(randperm(size(x,1),3000),:);
y=2.*mnorm_pdf(x, [0 0], [0.5 0;0 0.5]) + mnorm_pdf(x, [3 3], [0.5 0;0 0.5]) + mnorm_pdf(x, [-3 -3], [0.5 0;0 0.5]);
y=y+mnorm_pdf(x, [3 -3], [0.5 0;0 0.5])+mnorm_pdf(x, [-3 3], [0.5 0;0 0.5]);
y=y+0.03.*randn(size(y,1),1);
y(y>0.15)=1; y(y<=0.15)=-1;

[xt1, xt2]=meshgrid(-5:0.23:5,-5:0.23:5);
xt=[xt1(:) xt2(:)];
yt=ones(size(xt,1),1);

[n, nin] = size(x);

% Define covariance and likelihood functions and create model
gpcf = gpcf_sexp();
lik=lik_probit();
gp=gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-6);

opt=optimset('TolX',1e-3,'TolFun',1e-3,'display','on');
w0=gp_pak(gp);

% fPASS-GP with fixed size of 800 points in active set and data divided to
% 10 subsets with 4 sweeps over data
start=tic;[gp, indA]=passgp(gp, x, y, 'opt', opt, 'npass', 4, 'ninit', 800, 'nsub', 10, 'display', 'on', 'fixed', 'on', 'pexc', 0.1, 'optimn', 2);time=toc(start);
tt=time;
[Eft, Varft, lpyt, Eyt, Varyt]=gp_pred(gp, x(indA,:), y(indA,:), xt, 'yt', yt);
figure, [cc,hh]=contour(reshape(xt(:,1),size(xt1,1), size(xt1,1)), reshape(xt(:,2),size(xt1,1), size(xt1,1)), reshape(exp(lpyt),size(xt1,1), size(xt1,1)), [0.1 0.9]);
clabel(cc,hh); title('Pr(y==1) (fpass-gp)')
param=gp_pak(gp);

% PASS-GP with inclusion threshold 0.65, deletion threshold 0.99, intial
% size of 400 points in active set and 3 sweeps over data.
gp=gp_unpak(gp,w0);
start=tic;[gp, indA2]=passgp(gp, x, y, 'opt', opt, 'pinc', 0.65, 'pdel', 0.99, 'npass', 3, 'ninit', 400, 'nsub', 10, 'display', 'on', 'optimn', 2);time=toc(start);
tt2=time;
[Eft2, Varft2, lpyt2, Eyt2, Varyt2]=gp_pred(gp, x(indA2,:), y(indA2,:), xt, 'yt', yt);
figure, [cc,hh]=contour(reshape(xt(:,1),size(xt1,1), size(xt1,1)), reshape(xt(:,2),size(xt1,1), size(xt1,1)), reshape(exp(lpyt2),size(xt1,1), size(xt1,1)), [0.1 0.9]);
clabel(cc,hh); title('Pr(y==1) (pass-gp)')
param2=gp_pak(gp);

% Full Gaussian Process for comparison
gp=gp_unpak(gp,w0);
opt.Display='iter';
start=tic;gp=gp_optim(gp,x,y,'opt',opt);tt3=toc;
[Eft3, Varft3, lpyt3, Eyt3, Varyt3]=gp_pred(gp, x, y, xt, 'yt', yt);
figure, [cc,hh]=contour(reshape(xt(:,1),size(xt1,1), size(xt1,1)), reshape(xt(:,2),size(xt1,1), size(xt1,1)), reshape(exp(lpyt3),size(xt1,1), size(xt1,1)), [0.1 0.9]);
clabel(cc,hh); title('Pr(y==1) (full gp)')

% Display some statistics

mlpd_fpassgp=mean(mean(lpyt,2))
time_fpassgp=mean(tt)
mlpd_passgp=mean(mean(lpyt2,2))
time_passgp=mean(tt2)
mlpd_full=mean(lpyt3)
time_full=tt3

% Plot data and active sets for both methods
figure(4), subplot(1,2,1),  plot(x(y==1,1),x(y==1,2),'or',x(y==-1,1),x(y==-1,2),'ob'); 
hold all; plot(x(indA,1), x(indA,2), '.k'); title('Data and active set (fpass-gp)')
legend('y=1', 'y=-1', 'Active set for fpass-gp');
subplot(1,2,2),  plot(x(y==1,1),x(y==1,2),'or',x(y==-1,1),x(y==-1,2),'ob'); 
hold all; plot(x(indA2,1), x(indA2,2), '.k'); title('Data and active set (pass-gp)')
legend('y=1', 'y=-1', 'Active set for pass-gp');
setrandstream(prevstream);
