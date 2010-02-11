function [w, ix_eff, W, AX] = smlr_learning(label, X, Nfeat, varargin)
% Learning parameters of ARD-sparse multinomial logistic regression model.
%
% The estimation algorthm is derived from variational Bayesian method with 
% Laplace approximation in W-stap (SLR-Laplace).
% Note that the label vectors must be consisting of {0,1} value.
%
% -- Usage
% [w, ix_eff, W, AX] = smlr_learning(label, X, Nfeat, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [N*1]
% X     : Explanatory matrix                        [N*#feature]
% Nfeat : Number of features per class 
%
% -- Field of Optional Input
% w0    : Initial value of weight parameter w 
% ax0   : Initial value of relevance parameter ax
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
    help smlr_learning
end

%% # of parameters     
Nclass = length(unique(label));
Nparm  = Nfeat * Nclass;

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0'     , 'real'   ,  [],   ones(Nparm,1);...
     'nlearn'  , 'integer', [1 inf], 1000;...
     'nstep'   , 'integer', [1 inf], 100;...
     'amax'    , 'real'   , [0 inf], 1e8;...
     'wmaxiter', 'integer', [1 inf], 50;...
     'wdisplay', 'string' , {'iter', 'off', 'final', 'notify'}, 'final';...
     'gamma0'  , 'real'   , [0 inf], 0;...
     'isplot'  , 'boolean', []     , 0;...
     });

if ~isstruct(opt)
   error(opt);
end    
     
if length(opt.ax0) ~= Nparm
    error('The size of ax0 is incorrect !!');
end

% parameter 
ax0      = opt.ax0;
Nlearn   = opt.nlearn;
Nstep    = opt.nstep;
AMAX     = opt.amax;
wmaxiter = opt.wmaxiter;
wdisplay = opt.wdisplay;
gamma0   = opt.gamma0;
isplot   = opt.isplot; 

%---------------------
% VB 
%---------------------
w = zeros(Nparm,1);
ax = ax0;
ix_eff = [1:Nparm]';

% keep update
AX = opt.ax0;;
W = w;

for n = 1 : Nlearn
    
    ax_eff = ax(ix_eff);
    ax0_eff = ax0(ix_eff);
    w0_eff = w(ix_eff);

    [ixf, ixc] = ind2sub([Nfeat, Nclass], ix_eff);

    %%% W-step

    option = optimset('Gradobj','on','Hessian','on', 'Display', wdisplay, 'MaxIter', wmaxiter);
    [w_eff,f,eflag,output,grad,H]=fminunc(@linfunmlr, w0_eff, option,...
        label, ax_eff, X, ixf, ixc, Nclass);

    dS_eff = inv(H); % #parm*1
    %% A-step
%    ax_eff = (1-ax_eff.* diag(dS_eff))./(w_eff.^2);
        
    ax_eff = (1-ax_eff.* diag(dS_eff)+2*gamma0)./(w_eff.^2+2*gamma0./ax0_eff);

    %% Prune ineffective parameters
    w = zeros(Nparm,1);
    w(ix_eff) = w_eff;
    ax(ix_eff) = ax_eff;
    ix_eff = find(ax < AMAX);
    
    if mod(n, Nstep) == 0
        fprintf('Iterations : %d, Feature Remained: %d \n', n, length(ix_eff));
        AX = [AX ax];
        W  = [W  w];
        if isplot   % plot hyperparameters updating
            semilogy(ax+1);
            hold on
            semilogy(AMAX*ones(Nparm,1), 'r:', 'linewidth', 5)
            hold off
            title(['Iteration ', num2str(n), '  Feature ', num2str(length(ix_eff))]);
            pause(0.1);
        end
    end
end


