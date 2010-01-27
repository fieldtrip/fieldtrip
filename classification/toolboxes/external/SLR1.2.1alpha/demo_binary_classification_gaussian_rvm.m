% 2 class classfication by Gaussain RVM
% Nonlinear discriminant function can be used to classify low dimension features.
% Procedure which optimize the Gaussian kernel width using cross validation (CV) procedure .
% is shown. It should be noted a nonliear classifier may not work well in high dimensional
% problem. So it is required to reduce dimension before running this program.
%
% discriminant fucntion : kernel function 
% training data : simulated data generated from Gaussian Mixtures 
% test data     : simulated data generated from Gaussian Mixtures
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

clear
close all

fprintf('This is a demo how kernel classifier works with simulated data...\n\n');

%-----------------------
% Generate data
%-----------------------
mu = [-0.8 0.8; 0 0];
S  = [9 0; 0 1];
[tt, xx] = gen_sin_data(400,mu,S);
ttrain = tt(1:200);
xtrain = xx(1:200,:);
ttest  = tt(201:400);
xtest  = xx(201:400,:);

%------------------------
% Plot simulated data
%------------------------
slr_view_data(ttrain, xtrain);
xlabel('Feature 1');
ylabel('Feature 2');
title('Simulated Training Data');

fprintf('Press any key to proceed\n');
pause

%-------------------------
% Optimize Kernel Width
%-------------------------
fprintf('Kernel classifier with Gaussian basis function is used...\n'); 
fprintf('Now optimizing kernel width ...\n'); 

t = ttrain;
x = xtrain;
Nall = length(t);

RR = 4.^[1:10] * 1e-5; 
Ncv = 20;

for r = 1 : length(RR)
    R = RR(r);
    for nn = 1 : Ncv
        fprintf('Gaussian Width : %1.5f, Trial : %2.0d, ', R, nn);
        
        [ixtr, ixte] = separate_train_test(t, 0.7);
        t_train_CV = t(ixtr);
        t_test_CV  = t(ixte);
        x_train_CV = x(ixtr,:);
        x_test_CV  = x(ixte,:);
        
        [ww, ix_eff_all, errTable_tr, errTable_te, parm]...
            = biclsfy_rvm(x_train_CV, t_train_CV, x_test_CV, t_test_CV, R, ...
            'kernel', 'Gaussian', 'nlearn', 100, 'nstep', 150, 'scale_mode', 'none', 'mean_mode', 'none');

        COR_TEST(nn) = calc_percor(errTable_te);
          
    end  %% CV loop end

    COR_TEST_R(r,:) =  COR_TEST;
end

[tmp, rmax] = max(mean(COR_TEST_R, 2));

fprintf('Optimization finished ...\n'); 
fprintf('Optimum Kernel Width = %1.5f ...\n\n', RR(rmax)); 

%--------------------------
% Plot Optimization Result
%--------------------------

figure,
semilogx(RR, mean(COR_TEST_R, 2), 'bo-');
xyrefline(RR(rmax));   
xlabel('Kernel Width');
ylabel('Validation Result (%)');

fprintf('Press any key to proceed... \n');
pause

%-------------------------------
% Try Test Data 
%-------------------------------
fprintf('Now try test data ... \n');

Ropt = RR(rmax);

[ww, ix_eff_all, errTable_tr, errTable_te, parm]...
            = biclsfy_rvm(xtrain, ttrain, xtest, ttest, Ropt, ...
            'kernel', 'Gaussian', 'nlearn', 1000, 'nstep', 100, 'scale_mode', 'none', 'mean_mode', 'none');

            
fprintf('\nDemo Finished ...\n');
  
%---------------------------------
% Plot Boundary (need to modify)
%---------------------------------
xcenter = xtrain(ix_eff_all,:);    

minx1 = -10;
maxx1 = 10;
minx2 = -3;
maxx2 = 3 ;
    
[X1, X2] = meshgrid([minx1:(maxx1-minx1)/20:maxx1],[minx2:(maxx2-minx2)/20:maxx2]); 
Phi = slr_make_kernel([X1(:),X2(:)], 'Gaussian', xcenter, Ropt);

Z = Phi * ww(ix_eff_all);
        
%% 
ttest = label2num(ttest);
ix1 = find(ttest == 1);
ix2 = find(ttest == 2);

figure,
plot(xtest(ix1,1),xtest(ix1,2), 'b*');
hold on;
plot(xtest(ix2,1),xtest(ix2,2), 'r.');
hold on
contour(X1, X2, reshape(Z, [21,21]),[0,0]);
xlabel('Feature 1');
ylabel('Feature 2');
title('Test Data')
axis([-10 10 -3 3]);


        
