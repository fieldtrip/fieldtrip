% Example 6. Rigid CPD point-set registration. Full options intialization.
% if some options are not set, the defaults are used. 3D large stanford bunny point-set.
% Make use of Fast Gauss transform. Can be slow on older machines.
clear all; close all; clc;

load bun_zipper.mat; Y=X;

% add a random rigid transformation
R=cpd_R(rand(1),rand(1),rand(1));
X=rand(1)*X*R';

% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;          % show every iteration
opt.outliers=0.5;   % use 0.5 noise weight

opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default). Can be quite slow for large data sets.

opt.max_it=100;     % max number of iterations
opt.tol=1e-8;       % tolerance
opt.fgt=1;          % [0,1,2] if > 0, then use FGT. case 1: FGT with fixing sigma after it gets too small (faster, but the result can be rough)
                    %  case 2: FGT, followed by truncated Gaussian approximation (can be quite slow after switching to the truncated kernels, but more accurate than case 1)

 

% registering Y to X
Transform=cpd_register(X,Y,opt);

figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
