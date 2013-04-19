%DEMO_DERIVATIVEOBS  Regression problem demonstration with derivative 
%                    observations
%
%  Description
%    The regression problem consist of a data with one input variable,
%    two output variables with Gaussian noise; observations and 
%    derivative observations. The constructed model is full GP with
%    Gaussian likelihood.
%
%    The covariance matrix K includes now also covariances between
%    derivative observations and between derivative and latent
%    observations. With derivative observations, the K matrix is a
%    block matrix with following blocks:
%
%        K = [K_ll K_Dl'; K_Dl K_DD]
%
%    Where D refers to derivative and l to latent observation and
%       K_ll = k(x_i, x_j | th)
%       K_Dl = d k(x_i, x_j | th) / dx_i
%       K_DD = d^2 k(x_i, x_j | th) / dx_i dx_j
%
%    To include derivative observations in the inference:
%
%       - provide partial derivative observations in the
%       observation vector after output observations
%       y=[y;dy_1;...;dy_n]; for ex. if size(x)=[10 2] ->
%       size(y)=[30 1]
%
%       - gp_set(gp, 'derivobs', 'on')
%
%   The demo is organised in two parts:
%     1) data analysis without derivative observations
%     2) data analysis with derivative observations
%
%  See also  DEMO_REGRESSION1
%

% Copyright (c) 2010 Tuomas Nikoskinen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Create the data
tp=9;                                  %number of training points -1
x=[-2:4/tp:2]';
y=sin(x).*cos(x).^2;                   % The underlying process
dy=cos(x).^3 - 2*sin(x).^2.*cos(x);    % Derivative of the process
ns=0.06;                               % noise standard deviation

% Add noise
y=y + ns*randn(size(y));
% derivative observations are also noisy
dy=dy + ns*randn(size(dy));           
% observation vector with derivative observations
y2=[y;dy];

% test points
xt=[-3:0.05:3]';
nt=length(xt);

%========================================================
% PART 1 GP model without derivative obs
%========================================================
disp('GP model without derivative obs')

% Covariance function
pl = prior_t();
pm = prior_sqrtt();
gpcf = gpcf_sexp('lengthScale', 0.5, 'magnSigma2', .5, ...
                 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
% Use default Gaussian likelihood
gp = gp_set('cf', gpcf);

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'DerivativeCheck','on');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'opt',opt);
% Do the prediction
[Eft, Varft] = gp_pred(gp, x, y, xt);

% PLOT THE DATA

figure
%m=shadedErrorBar(p,Eft(1:size(xt)),2*sqrt(Varft(1:size(xt))),{'k','lineWidth',2});
subplot(2,1,1)
m=plot(xt,Eft,'k','lineWidth',2);
hold on
plot(xt,Eft+2*sqrt(Varft),'k--')
hold on
m95=plot(xt,Eft-2*sqrt(Varft),'k--');
hold on
hav=plot(x, y(1:length(x)), 'ro','markerSize',7,'MarkerFaceColor','r');
hold on
h=plot(xt,sin(xt).*cos(xt).^2,'b--','lineWidth',2);
%legend([m.mainLine m.patch h hav],'prediction','95%','f(x)','observations');
legend([m m95 h hav],'prediction','95%','f(x)','observations');
title('GP without derivative observations')
xlabel('input x')
ylabel('output y')

%========================================================
% PART 2 GP model with derivative obs
%========================================================
disp('GP model with derivative obs')

% Option derivobs set so that the derivatives are in use
gp = gp_set('cf', gpcf, 'derivobs', 'on');

% Set the options for the optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'DerivativeCheck','on');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y2,'opt',opt);
% Do the prediction
[Eft2, Varft2] = gp_pred(gp, x, y2, xt);
% Use predictions for function values only
Eft2=Eft2(1:nt);Varft2=Varft2(1:nt);

% PLOT THE DATA
% plot lines indicating the derivative

subplot(2,1,2)
m=plot(xt,Eft2,'k','lineWidth',2);
hold on
plot(xt,Eft2+2*sqrt(Varft2),'k--')
hold on
m95=plot(xt,Eft2-2*sqrt(Varft2),'k--');
hold on
hav=plot(x, y(1:length(x)), 'ro','markerSize',7,'MarkerFaceColor','r');
hold on
h=plot(xt,sin(xt).*cos(xt).^2,'b--','lineWidth',2);

xlabel('input x')
ylabel('output y')
title('GP with derivative observations')

i1=0;
a=0.1;
ddx=zeros(2*length(x),1);
ddy=zeros(2*length(x),1);
for i=1:length(x)
  i1=i1+1;
  ddx(i1)=x(i)-a;
  ddy(i1)=y(i)-a*dy(i);
  i1=i1+1;
  ddx(i1)=x(i)+a;
  ddy(i1)=y(i)+a*dy(i);
end

for i=1:2:length(ddx)
  hold on
  dhav=plot(ddx(i:i+1), ddy(i:i+1),'r','lineWidth',2);
end
legend([m m95 h hav dhav],'prediction','95%','f(x)','observations','der. obs.');

