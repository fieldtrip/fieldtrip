% Example 1. Affine CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear all; close all; clc;

load cpd_data2D_fish; Y=X;

% Add a random affine transformation
B=eye(2)+0.5*abs(randn(2,2));
X=X*B';


opt.method='affine';
Transform=cpd_register(X,Y,opt);

% Initial point-sets
figure,cpd_plot_iter(X, Y); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X, Transform.Y);  title('After');
