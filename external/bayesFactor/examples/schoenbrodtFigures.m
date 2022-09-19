function results = schoenbrodtFigures(nr,nrMC)
% Run the Bayes Factor Design Analysis simulations of Schoenbdort & Wagenmakers
%
% Schoenbrodt, F. D. & Wagenmakers, E. J. 
% Bayes factor design analysis: Planning for compelling evidence. 
% Psychon. Bull. Rev. 1–15 (2017). doi:10.3758/s13423-017-1230-y
% 
% Please note that the various numerical comparisons are not expected to
% match perfectly; there is randomization in the effects, in the noise
% added to the effects. Ballpark matches are all we can expect here and
% that seems to be the case. A direct comparison of bayesfFactor.ttest with
% the R package results does show an exact match.
%
% BK  -Jan 2019

if nargin <2
    nrMC = 20; % Quick - test only. Set to 10000 for reliable results
end
switch (nr)
    case 3
        %%
        % Figure 3 -  Design analysis for a two-sample T-test with N=20 and N=100
        % and a standardized effect size  of 0.5. From the text I infer
        % that they were looking at a one-tailed comparison. Because
        % designAnalysis uses  X<Y I pick 'left'.
        figure(2);
        clf;
        N = [20 100];
        results = bf.designAnalysis('N',N,'test','ttest2','sequential',false,'tail','left','effectSize',0.5,'scale',sqrt(2)/2,'nrMC',nrMC);  
        disp('***************')
        disp(['N = ' num2str(N)])
        disp(['False Positives (%): ' num2str(100*results.H0.pFalse)])% Under H0 above evidenceBoundary of 6
        disp(['False Negatives (%): ' num2str(100*results.H1.pFalse)]);% Under H1 below evidenceBoundary of 1/6
        disp(['True Positives (%) ' num2str(100*results.H1.pTrue)])% Under H1 aove evidenceBoundary of 6
        disp(['True Negatives (%) ' num2str(100*results.H0.pTrue)]);% Under H0 below evidenceBoundary of 1/6        
        
        %%
        % % Sample size determination to see how many samples are needed to
        % reach the evidence boundary in 95% (pSuccess) of cases under H1. (Rouder
        % finds 146 so for convenience we limit our search to a range
        % around that. In a real design you would explore a wider range.
        N = [100 120 140:1:150];
        pSuccess = 0.95;
        results = bf.designAnalysis('N',N,'test','ttest2','sequential',false,'tail','left','effectSize',0.5,'scale',sqrt(2)/2,'nrMC',nrMC,'pSuccess',pSuccess,'plot',false);  
        disp('**************')
        disp(['Necessary sample size for ' num2str(100*pSuccess) '% success under H1 is ' num2str(results.H1.N.min)]);       
        disp(['False Positives (%): ' num2str(100*results.H0.pFalse)])% Under H0 above evidenceBoundary of 6
        disp(['False Negatives (%): ' num2str(100*results.H1.pFalse)]);% Under H1 below evidenceBoundary of 1/6
        disp(['True negative  (%): '  num2str(100*results.H0.pTrue)]);
        
        
    case 4 
        %%
        % Figure 4 - sequential Bayes sampling. 
        % Also demonstrates the use of a distribution of a-prioir effect
        % sizes (instead ot a single predcited effect size)
        % Note that estimating the distribution of end points  would
        % require a larger number of MonteCarlo simulations (nrMC)- here we
        % just run 100 to illustrate the main point.
        effectFun =@(N) (0.5+0.1*randn([N 1])); % effect size is drawn from a distributionacross
        N = 20:1:200;  %     Sequential sampliung from 20 to 200 samples (truly open ended has not been implemented)
        figure(4);
        clf;
        results = bf.designAnalysis('N',N,'test','ttest2','sequential',true,'tail','left','effectSize',effectFun,'scale',sqrt(2)/2,'nrMC',nrMC,'plot',true);  
        disp(['False Positive Rate: ' num2str(100*results.H0.pFalse) '%'])
        disp(['False Negative Rate: ' num2str(100*results.H1.pFalse) '%'])
        disp(['Median Sample Size (H1): ' num2str(results.H1.N.median)])
        disp(['Median Sample Size (H0): ' num2str(results.H0.N.median)])
case 5
        %%
        % Figure 5 - Sequential sampling with asymmetric evidence
        % boundaries and a largger minimal N reduces false positives
        % Note that estimating the distribution of end points  would
        % require a larger number of MonteCarlo simulations (nrMC)- here we
        % just run 100 to illustrate the main point.
        
        effectFun =@(N) (0.5+0.1*randn([N 1])); % Effect sizes drawn from a normal distriburtion
        N = 40:1:100;  % Sample sizes 
        evidenceBoundary  = [1/6 30]; %Evidence for H0 is accepted at 1/6, evidence for H1 needs a BF of 30.
        figure(5);
        clf;
        results = bf.designAnalysis('N',N,'test','ttest2','sequential',true,'tail','left','effectSize',effectFun,'scale',sqrt(2)/2,'nrMC',nrMC,'plot',true,'evidenceBoundary',evidenceBoundary);  
        disp(['False Positive Rate: ' num2str(100*results.H0.pFalse) '%'])
        disp(['False Negative Rate: ' num2str(100*results.H1.pFalse) '%'])
        disp(['Median Sample Size (H1): ' num2str(results.H1.N.median)])
        disp(['Median Sample Size (H0): ' num2str(results.H0.N.median)])
end