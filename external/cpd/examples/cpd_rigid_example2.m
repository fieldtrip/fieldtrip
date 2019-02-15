% Example 2. Rigid CPD point-set registration. Full options intialization.
% if some options are not set, the defaults are used. 2D fish point-set.
clear all; close all; clc;

load cpd_data2D_fish; Y=X;

% add a random rotation
R=cpd_R(rand(1));
X=X*R';

% delete some points and add outliers.
X=[X(1:80,:); 0.5*randn(40,2)]-2;
Y=[Y(10:end,:); 0.5*randn(40,2)]-2;


% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;          % show every iteration
opt.outliers=0.6;   % use 0.6 noise weight
opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=1;      % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;     % max number of iterations
opt.tol=1e-8;       % tolerance


% registering Y to X
[Transform, Correspondence]=cpd_register(X,Y,opt);
% X(Correspondence,:) corresponds to Y

figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');

% registering X to Y
[Transform, Correspondence]=cpd_register(Y,X,opt);
cpd_plot_iter(Y, Transform.Y);  title('After registering X to Y');
