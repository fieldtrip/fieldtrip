function [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte] =...
    biclsfy_rlrlap(x_train, t_train, x_test, t_test, varargin)
% Binary classification by RLR with Laplace approximation (RLR-Lap).
%
% Normalization, parameter estimation, and performance evaluation are executed.
%
% -- Usage
% [ww, ix_eff_all, errTable_Tr, errTable_te, parm, AXall, Ptr, Pte] =...
%   run_smlr_reg(x_train, t_train, x_test, t_test,  varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode', 'string'  , {'all','each','none'}, 'all';...
%      'mean_mode' , 'string'  , {'all','each','none'}, 'all';...
%      'ax0'       , 'real'    , []                   ,  [];...
%      'nlearn'    , 'integer' , [1 inf]              ,  1000;...
%      'nstep'     , 'integer' , [1 inf]              ,  100;...
%      'amax'      , 'real'    , [0 inf]              ,  1e8;...
%      'wmaxiter'  , 'integer' , [1 inf]              ,  50;...
%      'wdisp_mode', 'string'  , {'iter', 'off', 'final', 'notify'}, 'iter';...
%      'isplot'    , 'boolean' , []                   ,  0;...
%      'usebias'   , 'boolean' , []                   ,  1;...
%      'norm_sep'  , 'boolean' , []                   ,  0;... 
%      'displaytext','boolean' ,  []                  ,  1;...
%     });
%
% --- Output
% ww          :   Estimated weight parameters. [Nfeat, Nclass]
% ix_eff_all  :   Index of features survived. cell array of {Nclass}
% errTable_tr :   Counting table of each label estimated. [Nclass, Nclass]
% errTbale_te :   Counting table of each label estimated. [Nclass, Nclass]
% parm        :   Parmaters used in this routine. [struct]
% AXall       :   History of hyperparameters updating. [Nfeat*Nclass Nlearn]
% Ptr         :   Probaility of observing every label in training data. [Nsamp_tr Nclass]
%                 This value is used to put a label on each sample.
% Pte         :   Probaility of observing every label in training data. [Nsamp_te Nclass]
%                 This value is used to put a label on each sample.
%
% 2009/08/10 OY
% * Bug fix when 't_test' has only a single label of 't_train'.
% 2006/10/23 OY
% * 'Nclass' is removed from inputs.
% 2006/09/12 OY
%  * A field "nprobe" is introduced.
% 2006/09/06 OY  
%  * Output format modified (error table as output)
%  * Comment modified
% 2006/08/02 OY
% 2006/05/26 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 4
    help biclsfy_rlrlap;
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
    {'scale_mode', 'string'  , {'all','each','stdall','stdeach','none'}, 'all';...
     'mean_mode' , 'string'  , {'all','each','none'}, 'all';...
     'ax0'       , 'real'    , [0 inf]              ,  1;...   % <--- scalar
     'nlearn'    , 'integer' , [1 inf]              ,  1000;...
     'nstep'     , 'integer' , [1 inf]              ,  100;...
     'amax'      , 'real'    , [0 inf]              ,  1e8;...
     'wmaxiter'  , 'integer' , [1 inf]              ,  50;...
     'wdisp_mode', 'string'  , {'iter', 'off', 'final', 'notify'}, 'iter';...
     'usebias'   , 'boolean' , []                   ,  1;...
     'norm_sep'  , 'boolean' , []                   ,  0;... 
     'displaytext','boolean' ,  []                  ,  1;...
     'gamma0'    , 'real'    , [0 inf]              ,  0;...
    });

if ~isstruct(parm)
   error(parm);
end

AMAX = parm.amax;
ax0 = parm.ax0; 
Nlearn = parm.nlearn;
Nstep = parm.nstep;
wdisp_mode = parm.wdisp_mode;
wmaxiter = parm.wmaxiter;
usebias  = parm.usebias;
norm_sep = parm.norm_sep;
displaytext   = parm.displaytext;
gamma0   = parm.gamma0;

%%
if displaytext
    fprintf('------------------------------------\n');
    fprintf('  Binary classification by RLR-LAP  \n');
    fprintf('------------------------------------\n');
end

% add a regressor for bias
if usebias == 1
    Nfeat = Nfeat+1;
end
Nparm = Nfeat;

% set ax0 
parm.ax0 = ax0;

% keep constant parameters
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
label = t_train - 1;

[ww, ix_eff, Wall, AXall] = rlr_learning(label, Xtr, @linfun,...
      'nlearn', Nlearn, 'nstep', Nstep,...
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

% 
errTable_tr = slr_error_table(t_train, t_train_est);
errTable_te = slr_error_table(t_test, t_test_est);

Pcorrect_tr = calc_percor(errTable_tr);
Pcorrect_te = calc_percor(errTable_te);

if displaytext,
fprintf(' Training Correct : %2.2f %%,  Test Correct : %2.2f %%\n', Pcorrect_tr, Pcorrect_te);
end
