% Make a point move in the 2D plane
% State = (x y xdot ydot). We only observe (x y).
% Generate data from this process, and try to learn the dynamics back.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)

ss = 4; % state size
os = 2; % observation size
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
H = [1 0 0 0; 0 1 0 0];
Q = 0.1*eye(ss);
R = 1*eye(os);
initx = [10 10 1 0]';
initV = 10*eye(ss);

seed = 1;
rand('state', seed);
randn('state', seed);
T = 100;
[x,y] = sample_lds(F, H, Q, R, initx, T);

% Initializing the params to sensible values is crucial.
% Here, we use the true values for everything except F and H,
% which we initialize randomly (bad idea!)
% Lack of identifiability means the learned params. are often far from the true ones.
% All that EM guarantees is that the likelihood will increase.
F1 = randn(ss,ss);
H1 = randn(os,ss);
Q1 = Q;
R1 = R;
initx1 = initx;
initV1 = initV;
max_iter = 10;
[F2, H2, Q2, R2, initx2, initV2, LL] =  learn_kalman(y, F1, H1, Q1, R1, initx1, initV1, max_iter);

