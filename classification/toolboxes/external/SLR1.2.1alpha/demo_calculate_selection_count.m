% Calculate SC-count value which can be used for feature selection
%
% discriminant fucntion : linear function 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

%
fprintf('\nThis is a demo how sparse logistic regression can be used to select meaningful features by cross validation...\n');

%----------- Parameters to be modified
Ncv = 100;     % # of cross validation
Rtrain = 0.8;  % Ratio of training data set 

%----------------------------
% Generate Data
%----------------------------

D = 300;
Ntr = 100; 
Nte = 100;
% mean
mu1 = zeros(D,1);
mu2 = [1.5; 0; zeros(D-2,1)];
% covariance
S = diag(ones(D,1));
ro = 0.8;
S(1,2) = ro;
S(2,1) = ro;


[ttr, xtr, tte, xte, g] = gen_simudata([mu1 mu2], S, Ntr, Nte);
              
fprintf('\nThe data is generated from 2 Gaussian Mixture model of which centers (mean) are different.\n'); 
fprintf('Input feature dimension is %d. \n',D);
fprintf('But only the first dimension has difference in mean value between two classes,\n'); 
fprintf('and the other dimension has same mean value.\n');
fprintf('Therefore only the first dimension is detected as a meaningful feature,\n') 
fprintf('if you select features by feature-wise t-value ranking method.\n');
fprintf('However due to the correlation between the second dimension and the first dimension,\n')
fprintf('inclusion of the second dimension makes classfication more accurate.\n');

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
slr_view_data(ttr, xtr);
axis equal;
title('Training Data')

fprintf('\nPress any key to proceed \n');
%pause

%----------------------------
% Cross Validation
%----------------------------
x = xtr;
t = ttr;

for nn = 1 : Ncv
    
    fprintf('\n\nCross Validation Trial : %3d \n', nn)
    
    [ix_train, ix_test] = separate_train_test(t, Rtrain);
            
    x_train = x(ix_train,:);
    t_train = t(ix_train);
    x_test = x(ix_test,:);
    t_test = t(ix_test);
    
    fprintf('\n\nFast version (ARD-Variational)!!\n')
    [ww, ix_eff, errTable_tr, errTable_te, g] = biclsfy_slrvar(x_train, t_train, x_test, t_test,...
        'nlearn', 200, 'nstep', 200, 'mean_mode', 'each', 'scale_mode', 'each');
        
    CVRes(nn).ix_eff_all = ix_eff;
    CVRes(nn).errTable_te = errTable_te;
    CVRes(nn).errTable_tr = errTable_tr;
    CVRes(nn).g           = g;    
end

% N-value 
SC = calc_SCval(CVRes, {'Survived'});
figure,
bar(SC);
ylabel('SC-value');
xlabel('Feature Index');



