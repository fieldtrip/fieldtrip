%  Nonrigid Example 6. Coherent Point Drift (CPD).
%  Nonrigid registration of 3D bunny point sets with use of FGT and Lowrank kernel approximation.

clear all; close all; clc;
load cpd_data3D_bunny.mat % load bunny set (8171x3) res with defformation

% Init full set of options %%%%%%%%%%
opt.method='nonrigid_lowrank'; % use nonrigid registration with lowrank kernel approximation
opt.numeig=30;                 % leave only 30 larges (out of 8171) eigenvectors/values to approximate G
opt.eigfgt=1;                  % use FGT to find the largest eigenvectore/values 

opt.beta=2;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=1;              % show every iteration
opt.outliers=0.1;       % use 0.7 noise weight
opt.fgt=2;              % use FGT to compute matrix-vector products (2 means to switch to truncated version at the end, see cpd_register)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-3;           % tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Transform, C]=cpd_register(X,Y, opt);


figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
