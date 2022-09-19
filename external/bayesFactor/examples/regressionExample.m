%% Example with continuous covariates (i.e. regression)

%% Liang et al.
% Compare the Liang et al formula (bf.bfFromR2) with the explciit
% integration implemented in the bf.anova function (so that it can handle
% continuous covariates as well as mixtures of continuous and categorical
% covariates).
nrSamples = 50;
% Create co-variates with different variance
nrRegressors = 2;
regressors = repmat([.1 3],[nrSamples 1]).*rand(nrSamples,nrRegressors);
% Add an interaction term
regressors = [regressors ,0.01*prod(regressors,2)];
noiseLevel = 1;    
nrIterations = 10; % Estimate N models with the same regressors, but different betas and noise.
bfLiang = nan(nrIterations,1);
bfFull = nan(nrIterations,1);
for i=1:nrIterations
    beta = randn(1,nrRegressors+1);
    y = regressors*beta'+noiseLevel*randn([nrSamples 1]);
    r = num2cell(regressors,1);
    T = table(y,r{:},'VariableNames',{'y','x1','x2','int'});    
    % Use the bf.regression function (which uses the Liang et al formula)
    bfLiang(i) = bf.regression(T,'y~x1*x2');
    % Use the more general bf.anova function
    bfFull(i)=bf.anova(T,'y~x1*x2','continuousScale',1);
end
% Show the direct comparison. This should be a slope 1 line.
figure(1);
clf;
subplot(1,2,1);
plot(bfLiang,bfFull,'.')
hold on
set(gca,'XScale','Log','YScale','Log')
plot(xlim,xlim,'k')
axis square;
xlabel 'BF Integrated'
ylabel 'BF Liang Formula'
subplot(1,2,2);
hist(bfFull./bfLiang)
xlabel 'BF Ratio (Anova/bfFromR2)'
ylabel '# Models'

%% Tooth Growth
% Analyze the tooth growth data set. It describes the tooth grwoth (len)
% after administering a dose of vitamin juice in the form of either orange
% juice or ascorbic acid (supp).
load toothgrowth
toothgrowth.doseAsFactor = categorical(toothgrowth.dose);
% Traditional analysis, using dose as a categorical variable.
tgCat = fitlm(toothgrowth,'len~supp*doseAsFactor');
anova(tgCat)
% Traditional analysis using log2(dose) as a continuous covariate.
toothgrowth.dose = log2(toothgrowth.dose);
tgCont = fitlm(toothgrowth,'len~supp*dose');
anova(tgCont)
% Now the corresponding Bayes Factors
bfFullCat = bf.anova(toothgrowth,'len~supp*doseAsFactor') % Full categorical model agains intercept
bfNoInteractionCat = bf.anova(toothgrowth,'len~supp+doseAsFactor'); % Without Interaciton 
bfFullCat/bfNoInteractionCat  % Bayes Factor that quantifies the evidence in favor of an interaction.
%
bfFullCont = bf.anova(toothgrowth,'len~supp*dose')  % Full continuous model against intercept-only,
bfNoInteractionCont = bf.anova(toothgrowth,'len~supp+dose'); % Without interaction.
bfFullCont/bfNoInteractionCont % Bayes Factor to quantify the evidence in favor of an interacion
% Compare the model that uses dose as a category with the model that uses dose as a continuous variable.
bfFullCat/bfFullCont  

%% Attitude
% Analyze the attitude data set  (See
% <https://richarddmorey.github.io/BayesFactor/#glm>)
% These are all continuous covariates, so could use bf.bfFromR2 too.
load attitude
% Compare the rating~complaints model with a set of alternatives. Here we
% just evaluate all the models relative to the intercept only alternative
% and then take the ratios. 
models = {'rating~complaints','rating~complaints+learning','rating~complaints+learning+advance','rating~complaints+raises','rating~complaints+privileges','rating~complaints+advance'};
bfAnova = bf.anova(attitude,models,'continuousScale',1);
% Show the comparisons
for i=1:numel(models)
    fprintf('%s:\t\t%3.2f\n',models{i},bfAnova(i)./bfAnova(1)); % Compare ~complaints with other models.
end
%%
% Alternatively (quicker), use the bf.regression function:
[bfRegression,lm] = bf.regression(attitude,models);
%The bfRegression and bfAnova calculations give the same result
pctError= mean(100*(bfRegression-bfAnova)./(bfRegression+bfAnova))


%% Test interactions
bfInt = bf.anova(attitude,'rating~complaints*raises','continuousScale',1)
bfIntReg= bf.regression(attitude,'rating~complaints*raises')




