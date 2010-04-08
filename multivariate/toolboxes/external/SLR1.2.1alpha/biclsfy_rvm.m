function [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte] =...
    biclsfy_rvm(x_train, t_train, x_test, t_test, R, varargin)
% Binary classification by Relevance Vector Machine (RVM).
% (slr_learning_var.m is used.)
% 
% -- Usage
% [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte] =
% biclsfy_rvm(x_train, t_train, x_test, t_test, R, varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
% R       : Width of Gaussian kernels. This value is not used when
% 'kernel=linear'.
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode', 'string' , {'all','each','none','stdall','stdeach'}, 'each';...
%      'mean_mode' , 'string' , {'all','each','none'}, 'each';...
%      'ax0'       , 'real'   ,  []                  ,  [];...
%      'nlearn'    , 'integer',  [1 inf]             ,  1000;...
%      'nstep'     , 'integer',  [1 inf]             ,  100;...
%      'amax'      , 'real'   ,  [0 inf]             ,  1e8;...
%      'usebias'   , 'boolean',  []                  ,  1;...
%      'norm_sep'  , 'boolean',  []                  ,  0;... 
%      'kernel'    , 'string' , {'linear', 'Gaussian'} , 'linear',...  
%      });
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
% 2009/05/29 OY
% * The function is renamed from 'run_rvm.m'.
%
% 2007/10/03 OY  Linear kernel is supported and used as a default kernel.
% 2006/11/16 OY  first version released
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 5
    help biclsfy_rvm;
    return
end

% char label -> number {1,2}
[t_train, label_names, Nclass] = label2num(t_train);
[t_test] = label2num(t_test, label_names);

[Nsamp_tr, Nfeat] = size(x_train);
Nsamp_te = size(x_test,1);



%% input check for optional parameter.
parm = finputcheck(varargin, ...
    {'scale_mode', 'string' , {'all','each','stdall','stdeach','none'}, 'each';...
     'mean_mode' , 'string' , {'all','each','none'}, 'each';...
     'ax0'       , 'real'   ,  []                  ,  [];...
     'nlearn'    , 'integer',  [1 inf]             ,  1000;...
     'nstep'     , 'integer',  [1 inf]             ,  100;...
     'amax'      , 'real'   ,  [0 inf]             ,  1e8;...
     'usebias'   , 'boolean',  []                  ,  1;...
     'norm_sep'  , 'boolean',  []                  ,  0;... 
     'kernel'    , 'string' , {'linear', 'Gaussian'}   , 'linear';...  
     'displaytext','boolean',  []                  ,  1;...
     });
 
if ~isstruct(parm)
   error(parm);
end

AMAX   = parm.amax;
ax0    = parm.ax0;
Nlearn = parm.nlearn;
Nstep  = parm.nstep;
usebias  = parm.usebias;
norm_sep = parm.norm_sep;
displaytext = parm.displaytext;
%

if displaytext
    fprintf('------------------------------------------\n');
    fprintf('   Binary classification by %s RVM   \n', parm.kernel);
    fprintf('------------------------------------------\n');
end

if usebias
    Nparm = Nsamp_tr+1;
else
    Nparm = Nsamp_tr;
end

% 
if isempty(ax0)
    ax0 = ones(Nparm,1);
    parm.ax0 = ax0;
end

parm.nclass = Nclass;
parm.nsamp_tr = Nsamp_tr;
parm.nsamp_te = Nsamp_te;
parm.nfeat    = Nfeat;
parm.nparm    = Nparm;
parm.R        = R;

% normalize (sacling and baseline addjustment)
if norm_sep == 0
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode, scale, base);
else
    [x_train, scale, base] = normalize_feature(x_train, parm.scale_mode, parm.mean_mode);
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode);
end

%----------------------------------
% RVM learning
%----------------------------------
if strcmp(parm.kernel, 'linear'),
    Phi_train = slr_make_kernel(x_train, parm.kernel);
elseif strcmp(parm.kernel, 'Gaussian'),
    Phi_train = slr_make_kernel(x_train, parm.kernel, [], R);
end

if usebias
    Xtr = [Phi_train, ones(Nsamp_tr,1)];
else
    Xtr = [Phi_train];
end

% Learning Step
label = t_train - 1 ; % {1,2} -> {0,1}
[ww, ix_eff, W, AX] = slr_learning_var2(label, Xtr, 'nlearn', Nlearn, 'nstep', Nstep, 'amax', AMAX);
AXall = AX;
ix_eff_all = ix_eff;

%----------------------
% Training Result
%----------------------
[t_train_est, Ptr] = calc_label(Xtr, [zeros(length(ww),1),ww]);

%---------------------
% Test Result
%---------------------    
if isempty(ix_eff)
    display('No kernel is survived ::: therefore no boundary in the feature space !');
    t_test_est = 2*ones(Nsamp_te,1);
    Pte = 1/2 * ones(Nsamp_te, 2);
else
    if usebias & ix_eff(end) > Nsamp_tr  % include bias term
        xcenter = x_train(ix_eff(1:end-1),:);
    else
        xcenter = x_train(ix_eff(1:end),:);   
    end

    if strcmp(parm.kernel, 'linear'),
        Phi_test = slr_make_kernel(x_test, parm.kernel, xcenter);
    elseif strcmp(parm.kernel, 'Gaussian'),
        Phi_test = slr_make_kernel(x_test, parm.kernel, xcenter, R);
    end

    if usebias & ix_eff(end) > Nsamp_tr  % include bias term
        Xte = [Phi_test ones(Nsamp_te,1)];
    else
        Xte = Phi_test;
    end

    [t_test_est, Pte] = calc_label(Xte, [zeros(length(ix_eff),1), ww(ix_eff)]);
end

errTable_tr = slr_error_table(t_train, t_train_est);
errTable_te = slr_error_table(t_test, t_test_est);

Pcorrect_tr = calc_percor(errTable_tr);
Pcorrect_te = calc_percor(errTable_te);

fprintf(' Training Correct : %2.2f %%,  Test Correct : %2.2f %%\n', Pcorrect_tr, Pcorrect_te);







