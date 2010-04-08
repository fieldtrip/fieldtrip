function [w, ix_eff, W, AX] = slr_learning_var2(label, X, varargin)
% Learning parameters of ARD-sparse logistic regression model (SLR-var).
%
% The estimation algorthm is derived from lower bound method. 
% The likelihood function is approximated by Gaussian distribution using variational parameters.  (SLR-var)  
%
% Note that label vectors must be consisting of {0,1} values. 
%
% -- Usage 
% [w, ix_eff, W, AX] = slr_learning_var2(label, X, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [#sample*1]
% X     : Explanatory matrix                        [#sample*#feature]
%
% -- Field of Optional Input
% ax0   : Initial value of relevance parameters ax
% xi0   : Initial value of variational parameters xi
% nlearn : # of learning 
% nstep  : # of step at which parameters updataing is kept.
% amax   : Truncation criteria. Parameters whose relevance paramater ax is l
%          arger than this value are eliminated from the further iterations.
% invhesssian : If 1, the inverse of Hessian (#feature*#feature) is used.
%               If 0, the diagonal elements of posterior variance is calculated without   
%               doing Hessian matrix inversion (Much faster when #feature >> #sample).   
% -- Output
% w_e    : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
% -- Example
% > [w_e, ix_eff, W, AX] = slr_learning_var(t, X,...
%   'nlearn', 100, 'nstep', 10);
%
% 
% 2009/05/30 OY
% * A field name of option, 'nohessian' is replaced with 'invhessian'.  
% 2007/12/28 By Okito Yamashita
% * Bug fix in theoritical calculation (line 120) --> must be check more
% carefully
%
% 2007/01/10 By Okito Yamashita
% * A field name 'fast' is changed to 'nohesssian', and also its variable type.
% 2006/10/04 by Okito Yamashita 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

%% Error check
if nargin < 2
    help slr_learning_var2
end

%% # of parameters     
[Nsamp,Nparm]  = size(X); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0',         'real',    [],      ones(Nparm,1);...
     'xi0',         'real',    [],      2*ones(Nsamp,1);...
     'nlearn',      'integer', [1 inf], 150;...
     'nstep',       'integer', [1 inf], 10;...
     'amax',        'real',    [0 inf], 1e8;...
     'invhessian',  'boolean', [],      1;...
     });
 
if ~isstruct(opt)
   error(opt);
end
    
if length(opt.ax0) ~= Nparm 
    error(['ax0 must be a vector of size ' num2str(Nparm), 'x 1 !!']);
end

if length(opt.xi0) ~= Nsamp 
    error(['xi0 must be a vector of size ' num2str(Nsamp), 'x 1 !!']);
end


%% Initial value for A-step and W-step
Nlearn = opt.nlearn;
Nstep  = opt.nstep;
AMAX   = opt.amax;
invhessian   = opt.invhessian;

W = [];
AX = [];

ax = opt.ax0;   
xi = opt.xi0;
w  = zeros(Nparm,1); 
ix_eff = [1:Nparm]'; %effective index 


for nlearning = 1 : Nlearn
    %% Effective parameters
    ax_eff = ax(ix_eff);
    w_eff = w(ix_eff);
    X_eff = X(:, ix_eff);
    
    %% W-step
    lam = tanh(xi/2)./xi/4;  

    
    if invhessian  % faster when Nsamp << Nfeat
        %H_eff = diag(ax_eff) + 2*X_eff'*diag(lam)*X_eff;
        H_eff = diag(ax_eff) + 2*X_eff'* (lam(:,ones(1,length(ix_eff))).*X_eff);
        S_eff = inv(H_eff);   % inverse of Nfeat*Nfeat
        % W-step
        w_eff = 1/2*S_eff*X_eff'*(2*label-1);
        
        %% A-step
        ax_eff = (1-ax_eff.*diag(S_eff))./(w_eff.^2);

        
        
      %  ax_eff = 1 ./ ((w_eff.^2)+diag(S_eff));
        %% Xi-step
        Xw = X_eff*w_eff;
        xi2 = Xw.^2 + sum(X_eff.*(X_eff*S_eff),2);
    else % faster when Nfeat >> Nsamp but slower when Nfeat == Nsamp
        ia = 1./ax_eff;
        ilam = 1./lam;
        iAX = ia(:,ones(1,Nsamp)) .* X_eff';
        XiAX = X_eff*iAX;

        C = 1/2*diag(ilam) + XiAX; % <-------------- 2007/12/28  OK 
        iC = inv(C);     % inverse of Nsamp*Nsamp 

        ib = zeros(length(ix_eff),1);
        for jj = 1 : length(ix_eff)
            ib(jj) = iAX(jj,:)*iC*iAX(jj,:)';
        end

        dS_eff = ia - ib;
        ilamy =   ilam.*(2*label-1);

        % W-step
        w_eff = 1/4* iAX * iC* ilamy;

        %% A-step
        ax_eff = (1-ax_eff.*dS_eff)./(w_eff.^2);
        %% Xi-step
        Xw = X_eff*w_eff;

        xi2 = Xw.^2 + 1/2*diag(XiAX*iC*diag(ilam));
        
    end
    %   %%%%%%%%%%%%%%%%%%%%
    %
    %   subplot(3,1,1),
    %   plot(dS_eff - diag(S_eff));
    %   subplot(3,1,2),
    %   plot(w_eff - w_eff2);
    %
    %   subplot(3,1,3)
    %   plot(xi2 - xi22);
    %     pause(2)
    
  
    
    xi = sqrt(xi2);
     
    %% Prune ineffective parameters
    w = zeros(Nparm,1);
    w(ix_eff) = w_eff;
    ax(ix_eff) = ax_eff;
    ix_eff = find(ax < AMAX);
  
  
   %% Keep history of parameter updating
   
   if mod(nlearning, Nstep) == 0
        fprintf('Iterations : %d, Feature Remained: %d \n', nlearning, length(ix_eff)); 
        W(:, nlearning/Nstep) = w;
        AX(:,nlearning/Nstep) = ax;
        
     %   semilogy(ax);
     %   pause(0.5);
   end
   
   if isempty(ix_eff)
       display('Caution: No feature is survived !!');
       break;
   end
   
end

AX = [opt.ax0 AX];
