function [w_eff, ix_eff, W, AX] = rlr_learning_var2(label, X, varargin)
% Learning parameters of ARD-regularized logistic regression model (RLR-var).
%
% The estimation algorthm is derived from lower bound method. 
% The likelihood function is approximated by Gaussian distribution using variational parameters.  (SLR-var)  
%
% Note that label vectors must be consisting of {0,1} values. 
%
% -- Usage 
% [w_eff, ix_eff, W, AX] = rlr_learning_var(label, X, varargin)
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
%
% -- Output
% w_e    : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
% AX     : History of hyper parameter learning
%
% -- Example
% > [w_e, ix_eff, W, AX] = rlr_learning_var2(t, X,...
%   'nlearn', 100, 'nstep', 10);
%
% 2009/05/29 OY
% * A field name of option is changed (nohessian --> invhessian).
% * Bug fix : inappropriate output arguments "w_eff" and "ix_eff"
% * Bug fix : ax_eff was not updated proprerly
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


%% Error check
if nargin < 2
    help rlr_learning_var2
end

%% # of parameters     
[Nsamp,Nparm]  = size(X); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'ax0',         'real',    [],      1;...
     'xi0',         'real',    [],      2*ones(Nsamp,1);...
     'nlearn',      'integer', [1 inf], 150;...
     'nstep',       'integer', [1 inf], 10;...
     'amax',        'real',    [0 inf], 1e8;...
     'invhessian',  'boolean', [],      1;...
     });
 
if ~isstruct(opt)
   error(opt);
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

ix_eff = [1:Nparm];
ax_eff = ax;
X_eff = X;

for nlearning = 1 : Nlearn
        
    %% W-step
    lam = tanh(xi/2)./xi/4;  

    if invhessian
        %H_eff = diag(ax_eff) + 2*X_eff'*diag(lam)*X_eff;
        H_eff = ax_eff*eye(Nparm) + 2*X_eff'* (lam(:,ones(1,Nparm)).*X_eff);
        S_eff = inv(H_eff);
        % W-step
        w_eff = 1/2*S_eff*X_eff'*(2*label-1);
        
        %% A-step
        %ax_eff = (Nparm-ax_eff.*diag(S_eff))./(w_eff.^2);
        ax_eff = Nparm ./ sum(w_eff.^2+diag(S_eff));
        %% Xi-step
        Xw = X_eff*w_eff;
        xi2 = Xw.^2 + sum(X_eff.*(X_eff*S_eff),2);
    else %%%%% New implementation
        ia = 1./ax_eff;
        ilam = 1./lam;
        iAX = ia*eye(Nparm)* X_eff';
        XiAX = X_eff*iAX;

        C = 1/2*diag(ilam) + XiAX;
        iC = inv(C);

        ib = zeros(length(ix_eff),1);
        for jj = 1 : length(ix_eff)
            ib(jj) = iAX(jj,:)*iC*iAX(jj,:)';
        end

        dS_eff = ia - ib;
        ilamy =   ilam.*(2*label-1);

        % W-step
        w_eff = 1/4* iAX * iC* ilamy;

        %% A-step
        %ax_eff = (Nparm-ax_eff.*diag(S_eff))./(w_eff.^2);
        ax_eff = Nparm ./ sum((w_eff.^2)+dS_eff);
       
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
     
    %% failure of regularization parameter estimation
   if ax_eff > AMAX
       display('Warning :: The regularization parameter can not be estimated correctly.');
       display('The result can not be reliable !!')
       w =zeros(Nparm,1),
       return
   end
     
   %% Keep history of parameter updating   
   if mod(nlearning, Nstep) == 0
        fprintf('Iterations : %d, Regularization Parameter: %f \n', nlearning, ax_eff); 
        W(:, nlearning/Nstep) = w_eff;
        AX(:,nlearning/Nstep) = ax_eff;
        
     %   semilogy(ax);
     %   pause(0.5);
   end
      
end

AX = [opt.ax0 AX];
