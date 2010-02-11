% Binary classfication demos
% discriminant fucntion : linear function 
%
% last updated 2009/06/15 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

data = 1;  % data = 1 or 2 or 3 (optional when you need to download testdata)

fprintf('This code demonstrates how a binary classification problem is solved ...\n');
%----------------------------
% Generate Data
%----------------------------
switch data
     case 1,
        D = 400;
        Ntr = 200;
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
        fprintf('But only the first dimension has difference in mean value between two classes,\n');
        fprintf('and the other dimension has same mean value.\n');
        fprintf('Therefore only the first dimension is detected as a meaningful feature,\n')
        fprintf('if you select features by feature-wise t-value ranking method.\n');
        fprintf('However due to the correlataion between the second dimension and the first dimension,\n')
        fprintf('inclusion of the second dimension makes classfication more accurate.\n');
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).\n');
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);

    case 2,
        D = 100;
        Ntr = 200;
        Nte = 100;
        % mean
        mu1 = zeros(D,1);
        mu2 = [[1:-0.02:0]'; zeros(D-51,1)];
        % covariance
        S = diag(ones(D,1));

        [ttr, xtr, tte, xte, g] = gen_simudata([mu1 mu2], S, Ntr, Nte);

        fprintf('\nThe data is generated from 2 Gaussian Mixture model of which centers (mean) are different.\n');
        fprintf('The mean value of the first 50 dimension is slightly different between two classes,\n');
        fprintf('whereas the remaining dimension has the same mean value.\n');
        fprintf('The degree of difference in the first 50 dimensions are manipulated\n')
        fprintf('by gradually changing the mean values in class 1 from 0 to 1.\n')
        fprintf('For comparison, this demo also computes the classification performance of linear RVM, \n');
        fprintf('which is Bayesian counterpart of support vector machine (SVM).\n');
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);
        
    case 3,
        load('../TESTDATA/real_binary', 'TRAIN_DATA', 'TEST_DATA', 'TRAIN_LABEL', 'TEST_LABEL');
        ttr = TRAIN_LABEL;
        tte = TEST_LABEL;
        xtr = TRAIN_DATA;
        xte = TEST_DATA;
        [Ntr,D] = size(TRAIN_DATA);
        [Nte] = size(TEST_DATA,1);
        
        fprintf('\nThe data is generated from a real experimental EEG data.\n');
        fprintf('In the experiment a subject executed either left or right finger tapping.\n');
        fprintf('The data has already been processed appropriately for classification.\n');
        fprintf('Input feature dimension is %d. \n', D); 
        fprintf('The number of training samples is %d. \n', Ntr);
end

%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------
slr_view_data(ttr, xtr);
axis equal;
title('Training Data')

fprintf('\n\nPress any key to proceed \n\n');
pause

%--------------------------------
% Learn Paramters
%--------------------------------
Algorithm = {'  SLR-LAP  ',...
             '  SLR-VAR  ',...
             '    RVM    ',...
             '  RLR-VAR  ',...
             'L1-SLR-LAP ',...
             'L1-SLR-COMP'};

tic
fprintf('\n%s!!\n', Algorithm{1});
[ww_o, ix_eff_o, errTable_tr(:,:,1) errTable_te(:,:,1)] = biclsfy_slrlap(xtr, ttr, xte, tte,...
    'wdisp_mode', 'off', 'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none');
time(1) = toc;
toc

tic
fprintf('\n%s!!\n', Algorithm{2});
[ww_f, ix_eff_f, errTable_tr(:,:,2), errTable_te(:,:,2)] = biclsfy_slrvar(xtr, ttr, xte, tte,...
    'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian',1);
time(2)=toc;
toc

tic
fprintf('\n%s!!\n', Algorithm{3});
    [ww_rvm, ix_eff_rvm, errTable_tr(:,:,3), errTable_te(:,:,3), g_rvm] = biclsfy_rvm(xtr, ttr, xte, tte, 0, ...
        'nlearn', 300, 'nstep', 100, 'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
time(3) = toc;
toc

tic
fprintf('\n%s!!\n', Algorithm{4});
    [ww_r, ix_eff_r, errTable_tr(:,:,4), errTable_te(:,:,4), g_r] = biclsfy_rlrvar(xtr, ttr, xte, tte, ...
        'nlearn', 300, 'nstep', 100, 'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
time(4) = toc;
toc

gamma = 10;
tic
fprintf('\n%s!!\n', Algorithm{5});  
[ww_l1, ix_eff_l1, errTable_tr(:,:,5), errTable_te(:,:,5), g_l1] = biclsfy_l1slrlap(xtr, ttr, xte, tte, sqrt(gamma),...
        'wdisp_mode', 'off', 'nlearn', 300, 'nstep', 100, 'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);
time(5) = toc;
toc

tic
fprintf('\n%s!!\n', Algorithm{6});
    [ww_l1c, ix_eff_l1c, errTable_tr(:,:,6), errTable_te(:,:,6), g_l1c] = biclsfy_l1slrc(xtr, ttr, xte, tte, gamma,...
        'nlearn', 300, 'nstep', 100, 'mean_mode', 'none', 'scale_mode', 'none');
time(6) = toc;
toc


% Weight vector conversion from kernels to features 
ww_rv = zeros(D,1);
for nn = 1 : Ntr,
ww_rv = ww_rv + ww_rvm(nn)*xtr(nn,:)';
end
ww_rv = [ww_rv; ww_rvm(Ntr+1)];

%--------------------------------
%  Result Table
%--------------------------------
fprintf('\n\n')
fmt1 = sprintf(' %%%ds       %%3.2f      %%3.2f      %%3.3f \n',11);
fmt2 = sprintf('                 %%5s    %%5s   %%6s\n');
fprintf(fmt2, 'Train(%)','Test(%)','Time(sec.)');
for ii = 1 : length(Algorithm),
    fprintf(fmt1,  Algorithm{ii}, calc_percor(errTable_tr(:,:,ii)), calc_percor(errTable_te(:,:,ii)), time(ii));
end
    



%--------------------------------
% Plot data (First 2 dimension)
%--------------------------------

figure,
subplot(3,2,1)
slr_view_data(tte, xte, [1 2], ww_o(:,1))
axis equal;
title('SLR-Laplce version');
subplot(3,2,2)
slr_view_data(tte, xte, [1 2], ww_f(:,1))
axis equal;
title('SLR-Variational version');
subplot(3,2,3)
slr_view_data(tte, xte, [1 2], ww_rv(:,1))
axis equal;
title('Linear RVM');
subplot(3,2,4)
slr_view_data(tte, xte, [1 2], ww_r(:,1))
axis equal;
title('RLR-Variational version');
subplot(3,2,5)
slr_view_data(tte, xte, [1 2], ww_l1(:,1))
axis equal;
title('L1-SLR-Laplace');
subplot(3,2,6)
slr_view_data(tte, xte, [1 2], ww_l1c(:,1))
axis equal;
title('L1-SLR-Component');



fprintf('\nFinish demo !\n');
