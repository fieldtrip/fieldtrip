function rouderFigure5
% Runs the simulations of Rouder et al. 2009 in Figure 5,
% to illustrate the differences between fixed and random effects.
load rouder2012Data
% Both ori and freq fixed
bfFullFixed= bf.anova(data,'rt~ori*freq');
bfBothFixed= bf.anova(data,'rt~ori+freq');
bfOriFixed= bf.anova(data,'rt~ori +ori:freq');
bfFreqFixed= bf.anova(data,'rt~freq+ ori:freq');

bf10  = nan(3,4); % 3 factors (ori,freq,int) and 4 kinds of effects (fixed, mixed, mixed, random)
bf10(1,1) = bfFullFixed/bfFreqFixed; % Main or int effect of ori
bf10(2,1) = bfFullFixed/bfOriFixed; % Main or int effect of freq
bf10(3,1) = bfFullFixed/bfBothFixed; % Interaction


% Ori fixed, freq random
bfFullMixed= bf.anova(data,'rt~ori*freq','treatAsRandom',{'freq'});
bfBothMixed= bf.anova(data,'rt~ori+freq','treatAsRandom',{'freq'});
bfOriMixed= bf.anova(data,'rt~ori +ori:freq','treatAsRandom',{'freq'});
bfFreqMixed= bf.anova(data,'rt~freq+ ori:freq','treatAsRandom',{'freq'});

bf10(1,2) = bfFullMixed/bfFreqMixed;
bf10(2,2) = bfFullMixed/bfOriMixed;
bf10(3,2) = bfFullMixed/bfBothMixed;


% Ori random , freq fixed
bfFullMixed= bf.anova(data,'rt~ori*freq','treatAsRandom',{'ori'});
bfBothMixed= bf.anova(data,'rt~ori+freq','treatAsRandom',{'ori'});
bfOriMixed= bf.anova(data,'rt~ori +ori:freq','treatAsRandom',{'ori'});
bfFreqMixed= bf.anova(data,'rt~freq+ ori:freq','treatAsRandom',{'ori'});

bf10(1,3) = bfFullMixed/bfFreqMixed;
bf10(2,3) = bfFullMixed/bfOriMixed;
bf10(3,3) = bfFullMixed/bfBothMixed;

% Both random
bfFullRandom =bf.anova(data,'rt~ori*freq','treatAsRandom',{'freq','ori'});
bfBothRandom= bf.anova(data,'rt~ori+freq','treatAsRandom',{'freq','ori'});
bfOriRandom= bf.anova(data,'rt~ori+ori:freq','treatAsRandom',{'freq','ori'});
bfFreqRandom= bf.anova(data,'rt~freq+ori:freq','treatAsRandom',{'freq','ori'});

bf10(1,4) = bfFullRandom/bfFreqRandom;
bf10(2,4) = bfFullRandom/bfOriRandom;
bf10(3,4) = bfFullRandom/bfBothRandom;


figure(5);
clf;

h = bar(1:3,bf10);
set(gca,'XTick',1:3,'XTIckLabel',{'Orientation','Frequency','Interaction'},'YScale','Log','YTick',[0.1 1 10 100 ])
set(gca,'Ylim',[0.1 200])
for i=1:4
    h(i).FaceColor=  (i-1)*0.3333*ones(1,3);
    h(i).BaseValue = 1;
end
legend(h,{'Both Fixed','Orientation Fixed','Frequency Fixed','Both Random'},'Location','NorthEast')
ylabel 'Bayes Factor vs. Null Model'

end