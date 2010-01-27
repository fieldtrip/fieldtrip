function [w, ix_eff, W, AX] = rmlr_learning(label, X, Nfeat, varargin)
% Learning parameters of regularized multinomial logistic regression model.
%
% The estimation algorthm is derived from variational Bayesian method with 
% Laplace approximation in W-stap (SLR-Laplace).
% Note that the label vectors must be consisting of {0,1} value.
%
% -- Usage
% [w, ix_eff, W, AX] = rmlr_learning(label, X, Nfeat, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [N*1]
% X     : Explanatory matrix                        [N*#feature]
% Nfeat : Number of features per class 
%
% -- Field of Optional Input
% w0    : Initial value of weight parameter w 
% ax0   : Initial value of relevance parameter ax
% prenlearn : # of prelearning 
% nlearn : # of learning 
% amax   : Truncation criteria. Parameters whose relevance parameter is larger 
%          than this value are eliminated from further iterations.
% wmaxiter : # of learning in W-step
% wdisplay : 
%            
% -- Output
% w      : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
% 2009/06/01 OY ver1.00
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


%% Error check
if nargin < 3
    help rmlr_learning
end

%% # of parameters     
Nclass = length(unique(label));
Nparm  = Nfeat * Nclass;

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0'      , 'real'   ,  [],   1;...
     'prenlearn', 'integer', [0 inf], 50 
     'nlearn'   , 'integer', [1 inf], 300;...
     'nstep'    , 'integer', [1 inf], 100;...
     'amax'     , 'real'   , [0 inf], 1e8;...
     'wmaxiter' , 'integer', [1 inf], 50;...
     'wdisplay' , 'string' , {'iter', 'off', 'final', 'notify'}, 'final';...
     'gamma0'   , 'real'   , [0 inf], 0;...
     });

if ~isstruct(opt)
   error(opt);
end    
     

% parameter 
ax0      = opt.ax0;
preNlearn = opt.prenlearn;
Nlearn   = opt.nlearn;
Nstep    = opt.nstep;
AMAX     = opt.amax;
wmaxiter = opt.wmaxiter;
wdisplay = opt.wdisplay;
gamma0   = opt.gamma0;

%---------------------
% VB 
%---------------------
w = zeros(Nparm,1);
axsc = ax0;    % <-- scalar

ix_eff = [1:Nparm]';
AX = ax0;
W = [];

for n = 1 : Nlearn
        
    ax_eff = axsc*ones(Nparm,1);  % vector
    w0_eff = w;

    [ixf, ixc] = ind2sub([Nfeat, Nclass], ix_eff);

    %%% W-step
    option = optimset('Gradobj','on','Hessian','on', 'Display', wdisplay, 'MaxIter', wmaxiter); 
    [w_eff,f,eflag,output,grad,H]=fminunc(@linfunmlr, w0_eff, option,...
             label, ax_eff, X, ixf, ixc, Nclass);

    dS_eff = inv(H); % #parm*1
    
    %% A-step
    if n < preNlearn  % additive rule (slower convergence)
        axsc = Nparm/(sum(diag(dS_eff))+sum((w_eff.^2)));
    else  % multiplicative rule (faster convergence)
        axsc = (Nparm - axsc * sum(diag(dS_eff)))/sum((w_eff.^2));
        %    ax_eff = (1-ax_eff.*
        %    diag(dS_eff)+2*gamma0)./(w_eff.^2+2*gamma0./ax0_eff);
    end

    if mod(n, Nstep) == 0,
       fprintf('Iterations : %d, Regularization Parameter: %f \n', n, axsc); 
       AX = [AX axsc];        
       W = [W w_eff];        
    end
    
   
    if axsc > AMAX
        w = zeros(Nparm, 1);
        break;
    end

    w = w_eff;
end




