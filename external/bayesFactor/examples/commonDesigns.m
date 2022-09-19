
% In a paper you'll often see a p-value and an N, based on either a T-test,
% a paired T-test or a one-way ANOVA.  To get an intuition how "good" a
% p-value is given the N, this script simulates some of those tests.

%% General parameter settings
nrSimulatedExperiments = 1000;
sampleSizes = [10 20 50];
testToSimulate = '2WAYRMMAIN';   % Pick one from 'TTEST','1WAY',2WAYRMMAIN 2WAYRMINT
nrSampleSizes= numel(sampleSizes);
scale =1;
bayesFactor = nan(nrSampleSizes,nrSimulatedExperiments);
bayesFactorBIC = nan(nrSampleSizes,nrSimulatedExperiments);
pValue = nan(nrSampleSizes,nrSimulatedExperiments);
stat =  nan(nrSampleSizes,nrSimulatedExperiments);
effectSize = 0.25;
tic;nn=0;

if ~exist('showTimeToCompletion.m','file')
    % Showtimetocompletion is a function outside the toolbox -do nothing
    showTimeToCompletion = @(x,y)(1); %NOP
end
switch upper(testToSimulate)
    case 'TTEST'
        %% Two-Sample T Test
        for i=1:nrSimulatedExperiments
            for j = 1:nrSampleSizes
                X = randn([sampleSizes(j) 1]);
                Y = randn([sampleSizes(j) 1])+effectSize;
                [bayesFactor(j,i),pValue(j,i),~,stats] = bf.ttest2(X,Y,'scale',scale);
                [bayesFactorBIC(j,i)] = bf.bfFromT(stats.tstat,stats.df); % Calculate BIC approximation
                stat(j,i) = stats.tstat;
                nn = showTimeToCompletion(((i-1)*nrSampleSizes+j)/(nrSimulatedExperiments*nrSampleSizes),nn);          
            end            
        end
    case '1WAY'
        %% 1-WAY ANOVA - main effect
        nrLevels = 4;
        cntr=0;
        for i=1:nrSimulatedExperiments
            for j = 1:nrSampleSizes
                X= repmat((0:nrLevels-1),[sampleSizes(j) 1]);
                X=X(:);
                y = randn([sampleSizes(j)*nrLevels 1])+X*effectSize;
                subject = repmat(1:sampleSizes(j),[nrLevels 1])';
                tbl = table(X,y,subject(:));
                [bayesFactor(j,i),lme] = bf.anova(tbl,'y~X','scale',scale);
                
                tmp =anova(lme);
                pValue(j,i) =tmp.pValue(2); % pValue of the Main factor
                stat(j,i) = tmp.FStat(2);
                df1 = nrLevels-1;
                df2 = nrLevels*(sampleSizes(j)-1);
                [bayesFactorBIC(j,i)] = bf.bfFromF(tmp.FStat(2),df1,df2,sampleSizes(j));
                nn = showTimeToCompletion(((i-1)*nrSampleSizes+j)/(nrSimulatedExperiments*nrSampleSizes),nn);                        end
        end
    case {'2WAYRMMAIN','2WAYRMINT'}
        %%
        levels = [2 3];
        nrFactors = numel(levels);
        for i=1:nrSimulatedExperiments
            for j = 1:nrSampleSizes
                fac1 = repmat((1:levels(1))',[1 levels(2)]);
                fac2 = repmat(1:levels(2),[levels(1) 1 ]);
                first = effectSize*randn([levels(1) 1]);
                second = effectSize*randn([1 levels(2)]);
                interaction = effectSize.*fac1.*fac2;
                singleSubject = repmat(first,[1 levels(2)]) + repmat(second,[levels(1) 1]) + interaction;
                response = reshape(repmat(singleSubject,[1 1 sampleSizes(j)])+ randn([levels sampleSizes(j)]),[prod(levels)*sampleSizes(j) 1]);
                fac1 = reshape(repmat(fac1,[1 1 sampleSizes(j)]),[prod(levels)*sampleSizes(j) 1]);
                fac2 = reshape(repmat(fac2,[1 1 sampleSizes(j)]),[prod(levels)*sampleSizes(j) 1]);
                subject = reshape(repmat((1:sampleSizes(j)),[prod(levels) 1 ]),[prod(levels)*sampleSizes(j) 1]);
                tbl = table(response,fac1,fac2,subject);
                switch upper(testToSimulate)
                    case '2WAYRMMAIN'
                        % Main effect Factor 1
                        full = bf.anova(tbl,'response~fac1+subject','scale',scale,'treatAsRandom',{'subject'});
                        restricted = bf.anova(tbl,'response~subject','scale',scale,'treatAsRandom',{'subject'});
                        bayesFactor(j,i) = full/restricted;
                        lm =fitlme(tbl,'response~fac1+(1|subject)');
                        tmp =anova(lm);
                        pValue(j,i) =tmp.pValue(2); % pValue of the Main factor
                        stat(j,i) = tmp.FStat(2);
                        df1 = levels(1)-1;
                        df2 = df1*(sampleSizes(j)-1);
                        bayesFactorBIC(j,i) = bf.bfFromF( stat(j,i),df1,df2,sampleSizes(j));
                        
                    case '2WAYRMINT'
                        % Interaction effect
                        full = bf.anova(tbl,'response~fac1*fac2+subject','scale',scale,'treatAsRandom',{'subject'});
                        restricted = bf.anova(tbl,'response~fac1+fac2+subject','scale',scale,'treatAsRandom',{'subject'});
                        bayesFactor(j,i) = full/restricted;
                        lm =fitlme(tbl,'response~fac1*fac2+(1|subject)');
                        tmp =anova(lm);
                        pValue(j,i) =tmp.pValue(end); % pValue of the Main factor
                        stat(j,i) = tmp.FStat(end);
                        df1 = prod(levels-1);
                        % df error for repeated measures: total-within-subjects
                        df2 = df1*(sampleSizes(j)-1);
                        bayesFactorBIC(j,i) = bf.bfFromF(stat(j,i),df1,df2,sampleSizes(j));
                end
                nn = showTimeToCompletion(((i-1)*nrSampleSizes+j)/(nrSimulatedExperiments*nrSampleSizes),nn);
            end
            %             consistency =100*mean(sign(log10(bayesFactor(j,:)))==sign(log10(bayesFactorBIC(j,:))));
            %             sprintf('%.2f %.2f %.2f %.2f %.2f %%',effectSize(1),min(log10(bayesFactor(:,j))),nanmedian(log10(bayesFactor(:,j))),max(log10(bayesFactor(:,j))),consistency)
            %             sprintf('%.2f %.2f %.2f %.2f %.2f %%',effectSize(1),min(log10(bayesFactorBIC(:,j))),nanmedian(log10(bayesFactorBIC(:,j))),max(log10(bayesFactorBIC(:,j))),consistency)
        end
        
end




%% Generate a figure
% 1 - compare pValue with BF
% 2 - compare stat (F or T) with BF
% 3 -  compare BIC approximation with BF

figure;clf
minP = 1e-4;
maxP = 0.1;
h = nan(nrSampleSizes,1);
xtick = sort([10.^(-3:1:-1) 0.05]);
for j = 1:nrSampleSizes
    subplot(2,2,1);
    x = pValue(j,:,:);
    [x,ix]=sort(x(:));
    y = bayesFactor(j,:,:);
    y=y(ix);
    h(j)= plot(x,y(:),'.-');
    set(gca,'YScale','Log','YLim',[0 1000],'XScale','log','XLim',[minP maxP],'XTIck',xtick);
    hold on
    subplot(2,2,2);
    x = stat(j,:,:);
    [x,ix]=sort(x(:));
    y = bayesFactor(j,:,:);
    y=y(ix);
    h(j)= plot(x,y(:),'.-');
    set(gca,'YScale','Log','YLim',[0 1000]);%,'XScale','log','XLim',[minP maxP],'XTIck',xtick);
    hold on
end
subplot(2,2,1)
plot(xlim,[1 1],'k--'); text(maxP,1,'Equal','Color','k')
plot(xlim,[3 3],'r--'); text(maxP,3,'Barely','Color','r')
plot(xlim,[6 6],'g--');text(maxP,6,'Moderate','Color','g')
plot(xlim,[10 10],'b--');text(maxP,10,'Strong','Color','b')
ylabel 'Bayes Factor'
xlabel 'p-Value'
title (['Simulating: ' testToSimulate])
bfLim =ylim;  % Take the bf limits that go with the p-limits we imposed
legend(h,num2str(sampleSizes'))
subplot(2,2,2);
plot(xlim,[1 1],'k--'); text(maxP,1,'Equal','Color','k')
plot(xlim,[3 3],'r--'); text(maxP,3,'Barely','Color','r')
plot(xlim,[6 6],'g--');text(maxP,6,'Moderate','Color','g')
plot(xlim,[10 10],'b--');text(maxP,10,'Strong','Color','b')
ylabel 'Bayes Factor'
xlabel 'Test Statistic (T/F)'
legend(h,num2str(sampleSizes'))

subplot(2,2,3)
for j = 1:nrSampleSizes
    x = bayesFactor(j,:,:);
    y = bayesFactorBIC(j,:,:);
    plot(x(:),y(:),'.');
    hold on
end
axis square
axis equal
ylabel 'BIC approximation'
xlabel 'Bayes Factor'
set(gca,'XScale','log','YScale','Log','XLim',bfLim,'YLim',bfLim);
plot(bfLim,bfLim,'k--')
title 'Comparing BF with BIC Approximation'