function [w, ix_eff, W] = slr_learning_l1c(label, X,  gamma1, varargin)
% Learning parameters of l1-sparse logistic regression model.
%
% The algorithm is based on Krishnapuram et.al 2005.
% Note that the label vectors must be consisting of {0,1} value.
%
% -- Usage
% [w, ix_eff, W] = slr_learning_l1c(label, X,  gamma1, varargin)
% 
% -- Input
% label : Teacher label vector consisting of {0,1}  [N*1]
% X     : Explanatory matrix                        [N*#feature]
% hfun  : Function handle
% gamma1 : a regularization parameter
%
% -- Field of Optional Input
% w0    : Initial value of weight parameter w 
% nlearn : # of learning 
% nstep  : # of steps to display updating results
%            
% -- Output
% w      : Estimated parameters
% ix_eff : Index of the effective feature
% W      : History of parameter learning
%
% -- Example
% >  [w, ix_eff] = slr_learning_l1c(t, X, gamma1, 'nlearn', 100, 'nstep', 10);
% The result of W-step is not displayed. # of iterations (W-step & A-step)
% is 100 and the updating hisotry of AX and W is kept in each 10 iterations. 
%
% 2009/06/10 OY 1st version
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

%% Error check
if nargin < 4
    help slr_learning_l1c
end

%% # of parameters     
Nparm  = size(X, 2); 

%% input check for optional parameter.
opt = finputcheck(varargin, ...
    {'w0',     'real',  [],   zeros(Nparm,1);...
     'nlearn', 'integer', [1 inf], 150;...
     'nstep',  'integer', [1 inf], 10;...
     });
    

if ~isstruct(opt)
   error(opt);
end    
     

%%
%% Initial Values
%%
w  = opt.w0; 
ix_eff = [1:Nparm]'; %effective index 

Nlearn = opt.nlearn;
Nstep  = opt.nstep;
Nsamp = length(label);

% Lower bound of Hessian 
Bl = -1/4*eye(Nsamp);
B = X'*Bl*X;

W = [];

for nlearning = 1 : Nlearn
  
   % Gradient
    y = X * w;
    p = 1 ./(1+exp(-y));
    gc = label-p;
    g = X'*gc;


   % component-wise update
   for dim = 1 : Nparm,
       w(dim) = soft(w(dim)-g(dim)/B(dim,dim), -gamma1/B(dim,dim)); 
   end
          
   ix_eff = find( w ~= 0);
%    if nlearning == 1, figure, end
%    semilogy(abs(w(ix_eff))); 
%    pause(0.1)
     
   
   %% Keep history of parameter updating
   
   if mod(nlearning, Nstep) == 0
       fprintf('Iterations : %d, Feature Remained: %d \n', nlearning, length(ix_eff)); 
       W(:, nlearning/Nstep) = w;      
   end
   
   if isempty(ix_eff)
       display('Caution: No feature is survived !!');
       break;
   end
   
end


%
% inner function
%

function x = soft(a,b)

c = abs(a) - b;

if c > 0
    x = c * sign(a);
else
    x = 0;
end
    
    



