function [w, ix_eff, W, AX] = slr_learning(label, X, hfun, varargin)
% Learning parameters of ARD-sparse logistic regression model.
%
% The estimation algorthm is derived from variational Bayesian method with 
% Laplace approximation in W-stap (SLR-Laplace).
% Note that the label vectors must be consisting of {0,1} value.
%
% -- Usage
% [w, ix_eff, W, AX] = slr_learning(label, X, hfun, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [N*1]
% X     : Explanatory matrix                        [N*#feature]
% hfun  : Function handle
%
% -- Field of Optional Input
% w0    : Initial value of weight parameter w 
% ax0   : Initial value of relevance parameter ax
% nlearn : # of learning 
% amax   : Truncation criteria. Parameters whose relevance parameter is larger 
%          than this value are eliminated from further iterations.
% wmaxiter : # of learning in W-step
% wdisplay : 
% reweight : If this is 'ON', the sparsity is accelarated in unusual way (no mathematical gurantee). 
%            
% -- Output
% w      : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
% -- Example
% >  [w, ix_eff] = slr_learning(t, X, @linfun,...
% > 'wdisplay', 'off', 'nlearn', 100, 'nstep', 10);
% The result of W-step is not displayed. # of iterations (W-step & A-step)
% is 100 and the updating hisotry of AX and W is kept in each 10 iterations. 
%
% 2006/10/24 by Okito Yamashita comment modified
% 2005/12/04 by Okito Yamashita 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

%% Error check
if nargin < 3
    help slr_learning
end

%% # of parameters     
Nparm  = size(X, 2); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0',    'real',  [],   ones(Nparm,1);...
     'w0',     'real',  [],   zeros(Nparm,1);...
     'nlearn', 'integer', [1 inf], 150;...
     'nstep',  'integer', [1 inf], 10;...
     'amax',   'real', [0 inf], 1e8;...
     'wmaxiter', 'integer', [1 inf], 50;...
     'wdisplay', 'string', {'iter', 'off', 'final', 'notify'}, 'final';...
     'reweight', 'string', {'ON', 'OFF'}, 'OFF'});

if ~isstruct(opt)
   error(opt);
end    
     
if length(opt.ax0) ~= Nparm | length(opt.w0) ~= Nparm
    error('The size of ax0 or w0 is incorrect !!');
end

%%
%% Initial value for A-step and W-step
%%
ax = opt.ax0;   
w  = opt.w0; 
ix_eff = [1:Nparm]'; %effective index 

Nlearn = opt.nlearn;
Nstep  = opt.nstep;
AMAX   = opt.amax;
ReWeight = opt.reweight;
WMaxIter = opt.wmaxiter;
WDisplay = opt.wdisplay;

W = [];
AX = [];

for nlearning = 1 : Nlearn
   %%
   %% Effective parameters 
   %%
   ax_eff = ax(ix_eff);
   w0_eff = w(ix_eff);
   X_eff = X(:, ix_eff);
   
   %%
   %% W-step 
   %%
   option = optimset('Gradobj','on','Hessian','on',...
       'MaxIter', WMaxIter, 'Display', WDisplay);
        
   [w_eff,f,eflag,output,g,H]=fminunc(hfun, w0_eff, option,...
       label, ax_eff, X_eff);
              
%    y = X_eff*w_eff;
%    p = 1 ./(1+exp(-y)) ; % #data
%    b = p.*(1-p);         % #data
%    
%    B = diag(b);       
%    A_eff = diag(ax_eff);
%    S_eff = inv(X_eff'*B*X_eff+A_eff);
 
   S_eff = inv(H);
  
   
   if strcmp(ReWeight, 'ON') 
   w_eff = S_eff*X_eff'*B* y;
   end

   %%
   %% A-step
   %%
   ax_eff = (1-ax_eff.*diag(S_eff))./(w_eff.^2);
   %ax_eff = 1./(w_eff.^2 + diag(S_eff));
 
   %%
   %% Prune ineffective parameters
   %%
   w = zeros(Nparm,1);
   w(ix_eff) = w_eff;
   ax(ix_eff) = ax_eff;
   ix_eff = find(ax < AMAX);
  
   
   %% Keep history of parameter updating
   
   if mod(nlearning, Nstep) == 0
       fprintf('Iterations : %d, Feature Remained: %d \n', nlearning, length(ix_eff)); 
       W(:, nlearning/Nstep) = w;
        AX(:,nlearning/Nstep) = ax;    
        
   end
   
   if isempty(ix_eff)
       display('Caution: No feature is survived !!');
       break;
   end
   
end

AX = [opt.ax0 AX];
