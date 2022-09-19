%%
% example from Richard Morey's Bayes Factor R-package
% http://bayesfactorpcl.r-forge.r-project.org/

sleep =   array2table( [1    0.7     1  1
2   -1.6     1  2
3   -0.2     1  3
4   -1.2     1  4
5   -0.1     1  5
6    3.4     1  6
7    3.7     1  7
8    0.8     1  8
9    0.0     1  9
10   2.0     1 10
11   1.9     2  1
12   0.8     2  2
13   1.1     2  3
14   0.1     2  4
15  -0.1     2  5
16   4.4     2  6
17   5.5     2  7
18   1.6     2  8
19   4.6     2  9
20   3.4     2 10],'VariableNames',{'Seq','extra','group','ID'});

% Analyze with Matlab mixed effects model  (Anova results in R are a
% bit different, not sure why).
lm = fitlme(sleep,'extra~group+(1+group:ID|ID)','FitMethod','REML');
a =anova(lm)
% Now do BF analysis
% First calculate BF for the model with group effect, and random effects ID
% and ID:group interaction (mimicking what R does with Error(ID/group)
[bfWithMain] = bf.anova(sleep,'extra~group+ID+group:ID','treatAsRandom',{'ID','group:ID'},'scale',1);
% Then fit a model with only the random effects
[bfRestricted] = bf.anova(sleep,'extra~ID+group:ID','treatAsRandom',{'ID','group:ID'},'scale',1);
% Take the ratio to quantify the evidence in favor of the main (fixed) effect of
% group.  Morey calculates the as ~11.6. I get 11.15 using gaussian
% quadrature. 
bfMainEffectOfGroup = bfWithMain/bfRestricted
