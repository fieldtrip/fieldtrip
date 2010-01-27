function [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte, pairs] =...
    muclsfy_slrvarovo(x_train, t_train, x_test, t_test, varargin)
% Multiclass classification by SLR-VAR one-versus-one classifier
% 
% -- Usage
% [ww, ix_eff_all, errTable_tr, errTable_te, parm, AXall, Ptr, Pte, pairs] =
% muclsfy_slrvarovo(x_train, t_train, x_test, t_test, varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode'  , 'string' , {'all','each','stdall','stdeach','none'}, 'each';...
%      'mean_mode'   , 'string' , {'all','each','none'}, 'each';...
%      'ax0'         , 'real'   ,  [],  [];...
%      'nlearn'      , 'integer',  [1 inf],  1000;...
%      'nstep'       , 'integer',  [1 inf],  100;...
%      'amax'        , 'real'   ,  [0 inf],  1e8;...
%      'usebias'     , 'boolean',  []     , 1;...
%      'norm_sep'    , 'boolean',  []     , 0;... 
%      'displaytext' , 'boolean',  []     , 1;... 
%      'invhessian'  , 'boolean',  []     , 0;...  
%      'combine_mode', 'integer', {0,1}   , 0;...
%      });
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
% pairs       :   Correspondence between linear index and pairs of binary classification  
%
%
% 2009/08/10 OY
% * Bug fix when 't_test' has only a single label of 't_train'.
% 2009/06/05 OY
% * the first version based on 'oy_Learn_SLR_pair'
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 4
    help  muclsfy_slrvarovo;
    return
end

% char label -> number 
[t_train, label_names, Nclass] = label2num(t_train);
[t_test] = label2num(t_test, label_names);

[Nsamp_tr, Nfeat] = size(x_train);
Nsamp_te = size(x_test,1);
Npair = Nclass*(Nclass-1)/2;

%% input check for optional parameter.
parm = finputcheck(varargin, ...
    {'scale_mode'  , 'string' , {'all','each','stdall','stdeach','none'}, 'each';...
     'mean_mode'   , 'string' , {'all','each','none'}, 'each';...
     'ax0'         , 'real'   ,  [],  [];...
     'nlearn'      , 'integer',  [1 inf],  1000;...
     'nstep'       , 'integer',  [1 inf],  100;...
     'amax'        , 'real'   ,  [0 inf],  1e8;...
     'usebias'     , 'boolean',  []     , 1;...
     'norm_sep'    , 'boolean',  []     , 0;... 
     'displaytext' , 'boolean',  []     , 1;... 
     'invhessian'  , 'boolean',  []     , 0;...  
     'combine_mode', 'integer', {0,1}   , 0;...
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
displaytext  = parm.displaytext;
invhessian   = parm.invhessian;
combine_mode = parm.combine_mode;

%
if displaytext
    fprintf('-----------------------------------------------------\n');
    fprintf('  Multi-class classification by SLR-VAR one-vs-one   \n');
    fprintf('-----------------------------------------------------\n');
end

% add bias
if usebias == 1
    Nfeat = Nfeat+1;
end

Nparm = Npair*Nfeat;


% initial parameters
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
    Xtr = [x_train];
    Xte = [x_test];
end

%---------------------
% one versus one 
%---------------------
ww = zeros(Nfeat,Npair);
ix_eff_all = cell(Npair,1);

AXall = [];
nn = 1;

for c1 = 1 : Nclass
    for c2 = c1+1 : Nclass
        if displaytext
            fprintf(': %d : Learning parameters for class %d vs class %d \n', nn, c1, c2);
        end

        ix1 = find( t_train == c1 );
        ix2 = find( t_train == c2 );

        Xpair_tr = Xtr([ix1(:);ix2(:)],:);
        label_pair = [zeros(length(ix1),1); ones(length(ix2),1)];  % 0 for c1, 1 for c2

        [ w_e, ix_eff, W, AX] = slr_learning_var2( label_pair, Xpair_tr, ...
            'nlearn', Nlearn, 'nstep', Nstep, 'amax', AMAX, 'ax0', ax0((nn-1)*Nfeat+1:nn*Nfeat),...
            'invhessian', invhessian);
        ix_eff_all{nn} = ix_eff;
        ww( :, nn ) = w_e;
        pairs(nn,:) = [c1 c2];
        AXall = [AXall; AX];
        nn = nn + 1;
    end
end

%-----------------------
% Training Correct
%-----------------------
[t_train_est, Ptr ] = calc_label_binary_pair( Xtr, ww, pairs, combine_mode );

%-----------------------
% Test
%----------------------
[t_test_est, Pte ] = calc_label_binary_pair( Xte, ww, pairs, combine_mode );

% remove bias parameters from effective indices
if usebias
    for cc = 1 : Npair
        ix_eff_all{cc} = setdiff(ix_eff_all{cc}, Nfeat);
    end
end

%
% Error Table (confusion matrix) and percent correct
%
errTable_tr = slr_error_table(t_train, t_train_est);
errTable_te = slr_error_table(t_test, t_test_est);

Pcorrect_tr = calc_percor(errTable_tr);
Pcorrect_te = calc_percor(errTable_te);

if displaytext,
fprintf(' Training Correct : %2.2f %%,  Test Correct : %2.2f %%\n', Pcorrect_tr, Pcorrect_te);
end