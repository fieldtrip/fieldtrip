%  Nonrigid Example 3. Coherent Point Drift (CPD).
%  Nonrigid registration of 3D face point sets with outliers.
%  Full set optioins is explicitelly defined. If you omit some options the
%  default values are used, see help cpd_register.

clear all; close all; clc;
load cpd_data3D_face.mat

% add outliers.
X=[X; randn(60,3)];
Y=[Y; randn(60,3)];

% Init full set of options %%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration
opt.beta=2;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=1;              % show every iteration
opt.outliers=0.4;       % noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-10;          % tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Transform, C]=cpd_register(X,Y, opt);


figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
