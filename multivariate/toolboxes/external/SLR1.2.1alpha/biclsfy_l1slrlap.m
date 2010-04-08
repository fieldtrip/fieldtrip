function [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall,Ptr,Pte] =...
    biclsfy_l1slrlap(x_train, t_train, x_test, t_test, gamma1, varargin)
% Binary classification by L1-norm-SLR with Laplace approximation (SLR-LAP).
% 
% -- Usage
% [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte] =
% biclsfy_slrlapl1(x_train, t_train, x_test, t_test, gamma1, varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
% gamma1  :   regualarization parameter (scalar) 
%             sqrt(gamma1) corresponds to the regularization parameter of
%             L1-norm
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode' ,'string' , {'all','each','stdall','stdeach','none'}, 'each';...
%      'mean_mode'  ,'string' , {'all','each','none'}, 'each';...
%      'ax0'        ,'real'   ,  []     ,  [];...
%      'nlearn'     ,'integer',  [1 inf],  1000;...
%      'nstep'      ,'integer',  [1 inf],  100;...
%      'amax'       ,'real'   ,  [0 inf],  1e8;...
%      'wmaxiter'   ,'integer',  [1 inf],  50;...
%      'wdisp_mode' ,'string' ,  {'iter', 'off', 'final', 'notify'}, 'iter';...
%      'usebias'    ,'boolean',  []     , 1;...
%      'norm_sep'   ,'boolean',  []     , 0;... 
%      'displaytext','boolean',  []     , 1;... 
%      });
%
% --- Output
% ww          :   Estimated weight parameters. [Nfeat, 1]
% ix_eff_all  :   Index of features survived. cell array of 1*1
% errTable_tr :   Counting table of each label estimated. [2, 2]
% errTbale_te :   Counting table of each label estimated. [2, 2]
% parm        :   Parmaters used in this routine. [struct]
% AXall       :   History of hyperparameters updating. [Nfeat Nlearn]
% Ptr         :   Probaility of observing every label in training data. [Nsamp_tr 2]
%                 This value is used to put a label on each sample.
% Pte         :   Probaility of observing every label in training data. [Nsamp_te 2]
%                 This value is used to put a label on each sample.
%
% 2009/08/10 OY
% * Bug fix when 't_test' has only a single label of 't_train'.
% 2009/06/09 OY  1st version
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 5
    help biclsfy_l1slrlap;
    return
end

% char label -> number 
[t_train, label_names, Nclass] = label2num(t_train);
[t_test] = label2num(t_test, label_names);

if Nclass ~= 2,
    error(' Use muclsfy_*.m !! ');
end

[Nsamp_tr, Nfeat] = size(x_train);
Nsamp_te = size(x_test,1);

%% input check for optional parameter.
parm = finputcheck(varargin, ...
    {'scale_mode' ,'string' , {'all','each','stdall','stdeach','none'}, 'each';...
     'mean_mode'  ,'string' , {'all','each','none'}, 'each';...
     'ax0'        ,'real'   ,  []     ,  [];...
     'nlearn'     ,'integer',  [1 inf],  1000;...
     'nstep'      ,'integer',  [1 inf],  100;...
     'amax'       ,'real'   ,  [0 inf],  1e8;...
     'wmaxiter'   ,'integer',  [1 inf],  50;...
     'wdisp_mode' ,'string' ,  {'iter', 'off', 'final', 'notify'}, 'iter';...
     'usebias'    ,'boolean',  []     , 1;...
     'norm_sep'   ,'boolean',  []     , 0;... 
     'displaytext','boolean',  []     , 1;... 
     });
 
if ~isstruct(parm)
   error(parm);
end
     
AMAX   = parm.amax;
ax0    = parm.ax0;
Nlearn = parm.nlearn;
Nstep  = parm.nstep;
wdisp_mode = parm.wdisp_mode;
wmaxiter = parm.wmaxiter;
usebias  = parm.usebias;
norm_sep = parm.norm_sep;
displaytext= parm.displaytext;

%
if displaytext
    fprintf('-------------------------------------- \n');
    fprintf('  Binary classification by L1-SLR-LAP  \n');
    fprintf('-------------------------------------- \n');
end
% add bias
if usebias == 1
    Nfeat = Nfeat+1;
end
% # of parameters
Nparm = Nfeat;

if isempty(ax0)
    ax0 = ones(Nparm,1);
    parm.ax0 = ax0;
end

parm.nclass = Nclass;
parm.nparm = Nparm;
parm.nsamp_tr = Nsamp_tr;
parm.nsamp_te = Nsamp_te;
parm.nfeat    = Nfeat;

% normalize (sacling and baseline addjustment)
if norm_sep == 0
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode, scale, base);
else
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode);
end

% add a regressor for bias term
if usebias
    Xtr = [x_train, ones(Nsamp_tr,1)];
    Xte = [x_test, ones(Nsamp_te,1)];
else
    Xtr = x_train;
    Xte = x_test;
end


%---------------------
% Learning
%---------------------
label = t_train - 1; % {1,2} --> {0,1}
         
[ww, ix_eff, W, AXall] = slr_learning_l1(label, Xtr, @linfun,gamma1,...
    'reweight', 'OFF', 'nlearn', Nlearn, 'nstep', Nstep,...
    'wdisplay', wdisp_mode, 'wmaxiter', wmaxiter, 'amax', AMAX);
     

%-----------------------
% Training Correct
%-----------------------
[t_train_est, Ptr] = calc_label(Xtr, [zeros(Nfeat,1) ww]);

%-----------------------
% Test
%----------------------
[t_test_est, Pte] = calc_label(Xte, [zeros(Nfeat,1) ww]);

% remove baseline parameters from effective index
if usebias
    ix_eff_all{1} = setdiff(ix_eff, Nfeat);
end

%---------------
% Performance
%---------------

errTable_tr = slr_error_table(t_train, t_train_est);
errTable_te = slr_error_table(t_test, t_test_est);

Pcorrect_tr = calc_percor(errTable_tr);
Pcorrect_te = calc_percor(errTable_te);

if displaytext,
fprintf(' Training Correct : %2.2f %%,  Test Correct : %2.2f %%\n', Pcorrect_tr, Pcorrect_te);
end