% KALMAN_FORWARD_BACKWARD		Forward Backward Propogation in Information Form 
%
%
% Note	:
%
%  M file accompanying my technical note
%
%    A Technique for Painless Derivation of Kalman Filtering Recursions
%
% available from http://www.mbfys.kun.nl/~cemgil/papers/painless-kalman.ps
% 

% Uses :

% Change History :
% Date		Time		Prog	Note
% 07-Jun-2001	 2:24 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

A = [1 1;0 1];
C = [1 0];
Q = eye(2)*0.01^2;
R = 0.001^2;
mu1 = [0;1];
P1 = 3*Q;

inv_Q = inv(Q);
inv_R = inv(R);

y = [0 1.1 2 2.95 3.78];

T = length(y);
L = size(Q,1);

%%%%% Forward message Passing 
h_f = zeros(L, T);
K_f = zeros(L, L, T);
g_f = zeros(1, T);
h_f_pre = zeros(L, T);
K_f_pre = zeros(L, L, T);
g_f_pre = zeros(1, T);


K_f_pre(:, :, 1) = inv(P1);
h_f_pre(:,1) = K_f_pre(:, :, 1)*mu1;
g_f_pre(1) = -0.5*log(det(2*pi*P1)) - 0.5*mu1'*inv(P1)*mu1;

for i=1:T,
  h_f(:,i) = h_f_pre(:,i) + C'*inv_R*y(:,i);
  K_f(:,:,i) = K_f_pre(:,:,i) + C'*inv_R*C;
  g_f(i) = g_f_pre(i) -0.5*log(det(2*pi*R)) - 0.5*y(:,i)'*inv_R*y(:,i);
  if i<T,
    M = inv(A'*inv_Q*A + K_f(:,:,i));
    h_f_pre(:,i+1) = inv_Q*A*M*h_f(:,i);
    K_f_pre(:,:,i+1) = inv_Q - inv_Q*A*M*A'*inv_Q;
    g_f_pre(i+1) = g_f(i) -0.5*log(det(2*pi*Q)) + 0.5*log(det(2*pi*M)) + 0.5*h_f(:,i)'*M*h_f(:,i);
  end;
end

%%% Backward Message Passing
h_b = zeros(L, T);
K_b = zeros(L, L, T);
g_b = zeros(1, T);

h_b_post = zeros(L, T);
K_b_post = zeros(L, L, T);
g_b_post = zeros(1, T);

for i=T:-1:1,
  h_b(:,i) = h_b_post(:,i) + C'*inv_R*y(:,i);
  K_b(:,:,i) = K_b_post(:,:,i) + C'*inv_R*C;
  g_b(i) = g_b_post(i) - 0.5*log(det(2*pi*R)) - 0.5*y(:,i)'*inv_R*y(:,i);
  if i>1,
    M = inv(inv_Q + K_b(:,:,i));
    h_b_post(:,i-1) = A'*inv(Q)*M*h_b(:,i);
    K_b_post(:,:,i-1) = A'*inv_Q*(Q - M)*inv_Q*A;
    g_b_post(i-1) = g_b(i) -0.5*log(det(2*pi*Q)) + 0.5*log(det(2*pi*M)) + 0.5*h_b(:,i)'*M*h_b(:,i);
  end;
end;


%%%% Smoothed Estimates

mu = zeros(size(h_f));
Sig = zeros(size(K_f));
g = zeros(size(g_f));
lalpha = zeros(size(g_f));

for i=1:T,
  Sig(:,:,i) = inv(K_b_post(:,:,i) + K_f(:,:,i));
  mu(:,i) = Sig(:,:,i)*(h_b_post(:,i) + h_f(:,i));
  g(i) = g_b_post(i) + g_f(:,i);
  lalpha(i) = g(i) + 0.5*log(det(2*pi*Sig(:,:,i))) + 0.5*mu(:,i)'*inv(Sig(:,:,i))*mu(:,i);
end;