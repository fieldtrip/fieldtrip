function rouderFigure2
% Runs the simulations of Rouder et al. 2009 in Figure 2,
% Shows the relationship between T test results, the Bayes Factor and
% sample size.
T = 1:0.05:10;
N = round(logspace(log10(5),log10(5000),100));
criticalT = nan(numel(N),3);
cntr =0;
for n=N
    cntr= cntr+1;
    for t= T
        % We call the bf.ttest function directly with
        % the results of a t-test (e.g. as returned by ttest.m)
        stats.tstat = t;
        stats.df    = n-1;
        stats.N     = n;  
        stats.tail  = 'both';
        stats.p     = NaN; % no tused
        % Call the class member
        bf10 = bf.ttest([],'stats',stats);
        % Check whether we;ve reach any of the thresholds.
        % If so, store.
        if bf10>3 && isnan(criticalT(cntr,1))
            criticalT(cntr,1) = t;
        end
         if bf10>10 && isnan(criticalT(cntr,2))
            criticalT(cntr,2) = t;
         end
         if bf10>30 && isnan(criticalT(cntr,3))
            criticalT(cntr,3) = t;
            break; % Goto next n
         end        
    end
end

%%
figure(2);
clf;
plot(N,criticalT)
xlabel 'Sample Size'
ylabel 'Critical t-value'
legend('BF=3','BF=10','BF=30')
set(gca,'XScale','Log','YLim',[1 6],'YTick',2:6,'XTick',[5 20 50 200 1000 5000],'XLim',[4 6000])
end
