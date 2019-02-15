% Example 1. Rigid CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear all; close all; clc;

load cpd_data2D_fish; Y=X;

% Add a random rotation and scaling
R=cpd_R(rand(1));
s=rand(1);
X=s*X*R';


Transform=cpd_register(X,Y);

% Initial point-sets
figure,cpd_plot_iter(X, Y); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X, Transform.Y);  title('After');

% Rotation and scaling errors after the registration
E_R=norm(R-Transform.R)
E_s=norm(s-Transform.s)
