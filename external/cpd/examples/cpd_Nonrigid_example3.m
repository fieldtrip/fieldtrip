%  Nonrigid Example 3. Coherent Point Drift (CPD).
%  Nonrigid registration of 3D face point sets with different deformations.
%  Full set optioins is explicitelly defined. If you omit some options the
%  default values are used, see help cpd_register.

clear all; close all; clc;
load cpd_data3D_face50.mat % load a model X and 50 deformed sets Y{}

% Init full set of options %%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration
opt.beta=2;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=1;              % show every iteration
opt.outliers=0;         % don't account for outliers
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-8;           % tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

face=50;                % try faces 1..50
[Transform, C]=cpd_register(X,Y{face}, opt);


X(:,1)=X(:,1)+3; % shift for side by side visualization
figure,cpd_plot_iter(X, Y{face}); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X.');
