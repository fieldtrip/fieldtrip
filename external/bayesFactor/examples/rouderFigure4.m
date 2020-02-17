function rouderFigure4(nrSets,options)
% Runs the simulations of Rouder et al. 2009 in Figure 4.
% This basically shows the ability to extract Main and Interaction effects
% in a 2-way anova.
if nargin<2
    options = bf.options;
    if nargin <1
        nrSets = 10;
    end
end
rouder2012 = load('rouder2012Data.mat');
effects = [0   0   0
    0.2 0   0
    0.5 0   0
    1   0   0
    0.2 0.4 0
    0.5 0.4 0
    1   0.4 0
    0.2 0.2   0
    0.5 0.5   0
    1   1     0
    0.4 0.4 0.2
    0.4 0.4 0.5];

nrEffects = size(effects,1);
% Create a design matrix from the data, only to simulate fake rt's with
% different effects
X= classreg.regr.modelutils.designmatrix(rouder2012.data,'intercept',false,'responsevar','rt','DummyVarCoding','effects','PredictorVars',{'ori','freq'},'model','interactions');
bfOri = nan(nrSets,nrEffects);
bfFreq= nan(nrSets,nrEffects);
bfInteraction= nan(nrSets,nrEffects);
parfor (j=1:nrEffects,options.nrWorkers)
    rt = X *effects(j,:)';
    for i=1:nrSets
        tmp  =rouder2012.data; %#ok<PFBNS>
        tmp.rt =rt + randn([size(X,1) 1]);
        bfFull = bf.anova(tmp,'rt~ori*freq');
        bfMain = bf.anova(tmp,'rt~ori+freq');
        bfOri(i,j)  = bf.anova(tmp,'rt~ori');
        bfFreq(i,j)  = bf.anova(tmp,'rt~freq');
        bfInteraction(i,j) = bfFull/bfMain;
    end
end

%%
figure(4);clf
effectIx = {1:4,5:7,8:10,11:12};
for i=1:4
    thisPlotEffects  = effectIx{i};
    subplot(1,4,i)
    if i==4
        x = effects(thisPlotEffects,3);
    else
        x = effects(thisPlotEffects,1);
    end
    [meOri]= median(bfOri(:,thisPlotEffects));
    hOri=  plot(x,meOri,'o-','MarkerSize',10,'Color','k','MarkerFaceColor','k');
    hold on
    meFreq = median(bfFreq(:,thisPlotEffects));
    hFreq=  plot(x,meFreq,'o-','MarkerSize',10,'Color','k','MarkerFaceColor',[0.7 0.5 0]);
    meInt = median(bfInteraction(:,thisPlotEffects));
    
    hInt=  plot(x,meInt,'o-','MarkerSize',10,'Color','k','MarkerFaceColor','w');
    if i==1
        legend([hOri, hFreq,hInt],{'Orientation','Frequency','Interaction'},'Location','NorthWest');
    end
    set(gca,'YScale','Log','YLim',[0.1 1e6],'YTick',[0.1 1 10 100 1000 1e6])
    ylabel 'Median Bayes Factor for Effect'
end
end