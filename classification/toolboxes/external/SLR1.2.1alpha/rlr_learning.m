function [w, ix_eff, W, AX] = rlr_learning(label, X, hfun, varargin)
% Learning parameters of regularized logistic regression model.
%
% The estimation algorthm is derived from variational Bayesian method with 
% Laplace approximation in W-stap (RLR-Laplace).
% Note that the label vectors must be consisting of {0,1} value.
%
% -- Usage
% [w, ix_eff, W, AX] = rlr_learning(label, X, hfun, varargin)
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
% 2009/05/29 OY
% * make a script this function
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


%% Error check
if nargin < 3
    help rlr_learning
end

%% # of parameters     
Nparm  = size(X, 2); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0',    'real',  [],   1;...
     'w0',     'real',  [],   zeros(Nparm,1);...
     'nlearn', 'integer', [1 inf], 150;...
     'nstep',  'integer', [1 inf], 50;...
     'amax',   'real', [0 inf], 1e8;...
     'wmaxiter', 'integer', [1 inf], 50;...
     'wdisplay', 'string', {'iter', 'off', 'final', 'notify'}, 'final';...
});

if ~isstruct(opt)
   error(opt);
end    
     

%
% Option 
%
Nlearn = opt.nlearn;
Nstep  = opt.nstep;
AMAX   = opt.amax;
WMaxIter = opt.wmaxiter;
WDisplay = opt.wdisplay;

%
% Initial value for A-step and W-step
%
ax0 = opt.ax0;   
w  = opt.w0; 

ax = ax0;    % <-- scalar
ix_eff = [1:Nparm];

W = [];
AX = [];

for n = 1 : Nlearn

    % variables to be updated :  ax, w
    
    ax_eff = ax*ones(Nparm,1);  % vector
    w0_eff = w;

    %%% W-step
    option = optimset('Gradobj','on','Hessian','on',...
        'Display', WDisplay, 'MaxIter', WMaxIter);
    
    [w_eff,f,eflag,output,grad,H]=fminunc(hfun, w0_eff, option,...
        label, ax_eff, X);
    dS_eff = inv(H); % #parm*1
    
    %% A-step
    %     if n > Nlearn/2
    %     axsc = (Nparm - axsc * sum(diag(dS_eff)))/sum((w_eff.^2));
    %     else
    %     axsc = Nparm/(sum(diag(dS_eff))+sum((w_eff.^2)));
    %     end

    %    ax_eff = (1-ax_eff.*
    %    diag(dS_eff)+2*gamma0)./(w_eff.^2+2*gamma0./ax0_eff);

    ax = Nparm/(sum(diag(dS_eff))+sum((w_eff.^2)));
   
    if mod(n, Nstep) == 0
       fprintf('Iterations : %d, Regularization Parameter: %f \n', n, ax); 
       AX(:, n/Nstep) = ax;        
       W(:, n/Nstep) = w_eff;        
    end

    
    if ax > AMAX
        display('Warning :: The regularization parameter can not be estimated correctly.');
        display('The result can not be reliable !!')
        w = zeros(Nparm, 1);
        break;
    end

    w = w_eff;   
end

AX = [ax0 AX];
