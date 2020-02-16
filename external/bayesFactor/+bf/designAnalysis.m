function [results] = designAnalysis(varargin)
% Perform a Bayes Factor Design Analysis.
%  This function will simulate data sets with a given effect
%  size, analyze  them using the test that is specified and
%  plot a historgram of BayesFactors under the null hypothesis
%  (effect size =0) and the H1 hypothesis (effect size is as
%  specified).
%
% Parm/Value pairs:
%
% N = Sample size: can be a vector
% sequential = Simulate a sequential Bayes Design (in this
%               design an increasing number of samples (N) is analyzed, and te
%               experiment stops when the Bayes Factor reaches a
%               pre-specified thredhold for either H1 or H0.
%               [false].
% evidenceBoundary = Boundary at which to stop collecting more
%                   samples. This should contain two elements,
%                   the threshold below which the evidence for H0 is
%                   considered sufficient and the threshold above
%                  which the evidence for H1 is considered sufficient. [1/6 6].
% nrMC  =   The number of monte carlo simulations to run [10000]
% test  = Which statistical test to use. TTEST, TTEST2,
% effectSize = The (expected) effectSize under H1. The "units" differ per test:
%                   TTEST,TTEST2: expressed as fraction of the standard deviation. [0.5]
%                   To specify a distribution of effectSizes,
%                   use an anonymouys function. For instance,
%                   for an effect size with a mean of 0.5 and stdev of 0.1,
%                   pass @(N) (0.5+0.1*randn([N 1]])
%
% tail = Tail for TTest/TTest2. [both]
% scale = Scale of the Cauchy prior [sqrt(2)/2].
% pSuccess = The target probability of "success" (reaching the
% evidence boundary). [0.9]
% plot = Toggle to show graphical output as in Schoenbrodt &
%                   Wagenmakers. Pass the results struct of a previous run in this argument
%                   to generate only the plot.
% options = Monte Carlo /Parralel execution options [bf.options]
%
% Using Linear Mixed Models
% linearModel = The linear model to simulate.
% alternativeModel = A Wilcoxon formula describing the alternative model to
% test against. [Defaults to the intercept only model See bf.anova)]
% subject = The variable in the model that identifies unique subjects. (see
% bf.internal.simulateLinearModel)
% In a linear model simulation, effectSize is the scaling factor of the
% standard deviation of the residuals in the simulated data.
%
% To add your own analysis, but reuse the looping/graphing from
% this funcition, set 'test' to 'SPECIAL' and provide
% 'dataFun', functions that generate data
% based on an effect size and N, and 'bfFun', a function that
% cacluates the bayesFactor based on the output of the dataFun.
% In code the following lines will be executed with the functions you
% provide
%  data  = dataFun(effectSize,N);
%  thisBf  =  bfFun(data{:});
%
% For sequential designs, there are some additional hurdles for dataFun ouput (see
% LinearMixedModel in the code below).
%
% For a tutorial on BFDA, see
%           Schoenbrodt, F. D. & Wagenmakers, E. J.
%           Bayes factor design analysis: Planning for compelling evidence.
%           Psychon. Bull. Rev. 1?15 (2017). doi:10.3758/s13423-017-1230-y
%
% OUTPUT
% A struct contining structs H1 and H0. They contain information on the simulations
% under H1 and H0, respectively. The content is different when
% a sequential verus a fixed design is used:
% FIXED N Design
%   .bf.all = The BayesFactor for each simulation and each N [nrMC nrDifferentSampleSizes]
%   .bf.median = Median bayesfactor, for each N.
%   .N.min  = This is the smalest sample size (out of the vector provided as input 'N') for which the
%               probability of the evidence reaching the evidence boundary is 'pSuccess' (an input parm)
%   .N.median = NaN (not used in Fixed-N)
%   .N.pMax = NaN (not used in Fixed-N)
% .pFalse - Propbabilty of misleading evidence (false positive
% for H0, false negative for H1).
% SEQUENTIAL design
%  .bf.all = The bayes factor for each simulation and each
%  sequential sample size . Each row represents an experiment
%  with increasing sample sizes. The sequence stops when an
%  evidence boundary is reached.  (bf will be NaN for larger
%  sample sizes in that row.)
% .bf.median = NaN = not used in sequential
% .N.all = The distribution of sample sizes
% N.median  = the mean sample size used in the sequential
%               design (median across all sequences that hit an
%               evidence boundary an dthose that reached n-max)
% N.pMax = probablity of terminatinng at n-max (i.e. not
% hitting an evidence boundary before N =n-Max).
% N.min = NaN - not used in a sequential design.            %
% % .pFalse - Propbabilty of misleading evidence (false positive
% for H0, false negative for H1).
%
% EXAMPLES:
%  Explore the power of a one sample T test with 20 or 100 subjects:
%   This will create a figure like Figure #3 in Schoenbrodt &
%   Wagenmakers.
%       bf.designAnalysis('N',[20 100],'test','ttest','sequential',false);
%
% Investigate a sequential sampling design where we obtain
% up to 200 samples but stop data collection if a BF10 of 1/6
% or 6 is reached. (The example illustrates that sample size
% spacing does not have to be regular
%    bf.designAnalysis('N',[20:2:60 65:5:120 130:10:200],'evidenceBoundary',[1/6 6],'test','ttest','sequential',true,'nrMC',1000);
%
% BK - Jan 2019

p = inputParser;
p.addParameter('N',[20 100]); % A vector of sample sizes to evaluate
p.addParameter('sequential',false); % Toggle to use a sequential design [false]
p.addParameter('nrMC',10000); % Each N is evaluated this many time (Monte Carlo)
p.addParameter('evidenceBoundary',[1/6 6]);
p.addParameter('pSuccess',0.9); % For a fixed-N design, which fraction of success (i.e. reaching the h1 evidence boundary) is acceptable (used to calculate minimalN)
p.addParameter('test','TTEST');  % Select one of TTEST, TTEST2,SPECIAL
p.addParameter('effectSize',0.5);
p.addParameter('linearModel',[]); % Used only for linearModel
p.addParameter('alternativeModel',''); % Used only for linearModel
p.addParameter('tail','both');
p.addParameter('scale',sqrt(2)/2);
p.addParameter('plot',true);
p.addParameter('bfFun',[]);  % For 'test'=='SPECIAL' - specify your own bfFun and dataFun here.
p.addParameter('dataFun',[]);
p.addParameter('options',bf.options);
p.parse(varargin{:});
results.H1 = struct('N',struct('min',NaN,'median',NaN,'all',[],'pMax',NaN),'bf',struct('all',[],'median',nan),'pFalse',[]);
results.H0 = results.H1;

if ~isempty(p.Results.linearModel)
    test = 'LINEARMODEL';
else
    test = p.Results.test;
end
% Depending on the test, create anonymous functions that will
% generate simulated data aand call the appropriate BF function
% in the MonteCarlo loop below.
switch upper(test)
    case 'TTEST'
        % BFDA for a one-sample ttest.
        dataFun = @(effectSize,N) ({effectSize + randn([N 1])});
        bfFun = @(X) bf.ttest(X,'tail',p.Results.tail,'scale',p.Results.scale);
    case 'TTEST2'
        % Two sample ttest
        dataFun = @(effectSize,N) ({-0.5*effectSize + randn([N 1]),0.5*effectSize + randn([N 1])});
        bfFun = @(X,Y) bf.ttest2(X,Y,'tail',p.Results.tail,'scale',p.Results.scale);
    case 'LINEARMODEL'
        %Linear mixed model analysis (using bf.anova)
        dataFun = @(effectSize,N) bf.internal.simulateLinearModel(p.Results.linearModel,effectSize,N);
        bfFun   = @(X,Y) bf.anova(X,Y,'alternativeModel',p.Results.alternativeModel);
    case 'SPECIAL'
        % User specified functions to generate simulated data
        % and calculate the bf
        dataFun = p.Results.dataFun;
        bfFun = p.Results.bfFun;
    otherwise
        errror(['BFDA for ' p.Results.testFun ' has not been implemented yet'])
end

% Initialize
N = p.Results.N;
nrN = numel(N);
results.H1.bf.all = nan(p.Results.nrMC,nrN);
results.H0.bf.all = nan(p.Results.nrMC,nrN);
maxN = max(p.Results.N);
upperBF = max(p.Results.evidenceBoundary);
lowerBF = min(p.Results.evidenceBoundary);

if isa(p.Results.effectSize,'function_handle')
    es = p.Results.effectSize(maxN);
else
    es = p.Results.effectSize;
end

%% Do the actual computation
if isstruct(p.Results.plot)
    showGraph = true;
    % Plot a previously computed results struct
    results =p.Results.plot;
else
    showGraph = p.Results.plot;
    if p.Results.sequential
        % Run a sequential design :  generate fake data for maxN
        % subjects but sample htese successively and at each step
        % test whether the BF is above threshold (and if so,
        % terminate). Repeat this nrMC times.
        nrFalsePositive = 0;
        nrFalseNegative = 0;
        h1BfAll = nan(p.Results.nrMC,nrN);
        h0BfAll = nan(p.Results.nrMC,nrN);
        parfor (j=1:p.Results.nrMC,p.Results.options.nrWorkers)
        %for j=1:p.Results.nrMC
            for hyp=0:1
                thisBf = nan(1,nrN);
                if hyp==0
                    % H0
                    data = dataFun(0,maxN);
                else
                    % H1
                    data  = dataFun(es,maxN);
                end
                for i=1:nrN
                    % Sequential sampling for H1 data
                    if strcmpi(test,'LinearModel')
                        %data{1} is the table. , data{2} the formula.
                        perN = height(data{1})/maxN; % For N=1, pick perN from the table.
                        thisData = {data{1}(1:perN*N(i),:),data{2}};
                    else
                        thisData = cellfun(@(x)(x(1:N(i))),data,'UniformOutput',false);
                    end
                    thisBf(i) = bfFun(thisData{:});
                    if thisBf(i) > upperBF || thisBf(i) < lowerBF
                        % Crossed threshold. Stop sampling.
                        if thisBf(i) <lowerBF && hyp==1
                            nrFalseNegative= nrFalseNegative+1;
                        else  % thisBf>upperBf
                            if hyp==0
                                nrFalsePositive= nrFalsePositive+1;                            
                            end
                        end
                        break;
                    end
                end
                if hyp==0
                    h0BfAll(j,:) = thisBf;
                else
                    h1BfAll(j,:) = thisBf;
                end
            end
        end
        results.H1.bf.all= h1BfAll; % Assign outside parfor
        results.H0.bf.all = h0BfAll;
        
        % Analyze the results to determine how many samples will
        % be collected
        [ix,col] = find(isnan(results.H1.bf.all)');
        [~,isFirst] =unique(col);
        nrSamplesEarlyStop  = p.Results.N(ix(isFirst));
        nrMax = p.Results.nrMC- numel(isFirst);
        results.H1.N.all = [nrSamplesEarlyStop maxN*ones(1,nrMax)];
        results.H1.N.median = median(results.H1.N.all);
        results.H1.pFalse  = nrFalseNegative./p.Results.nrMC;
        results.H1.N.pMax = nrMax./p.Results.nrMC;
        
        [ix,col] = find(isnan(results.H0.bf.all)');
        [~,isFirst] =unique(col);
        nrSamplesEarlyStop  = p.Results.N(ix(isFirst));
        nrMax = p.Results.nrMC- numel(isFirst);
        results.H0.N.all = [nrSamplesEarlyStop maxN*ones(1,nrMax)];
        results.H0.N.median = median(results.H0.N.all);
        results.H0.pFalse  = nrFalsePositive./p.Results.nrMC;
        results.H0.N.pMax = nrMax./p.Results.nrMC;
        
        results.H0.pTrue =  (p.Results.nrMC-nrMax-nrFalsePositive)/p.Results.nrMC;
        results.H1.pTrue =  (p.Results.nrMC-nrMax-nrFalseNegative)/p.Results.nrMC;
        
    else
        % Simualte a regular fixed N design
        h1BfAll = nan(p.Results.nrMC,nrN);
        h0BfAll = nan(p.Results.nrMC,nrN);
        parfor (j= 1:p.Results.nrMC,p.Results.options.nrWorkers)
            thisBf0 = nan(1,nrN);
            thisBf1 = nan(1,nrN);
            for i =1:nrN
                data  = dataFun(es,N(i));
                thisBf1(i) =  bfFun(data{:}); % Collect BF under H1
                nullData  = dataFun(0,N(i));
                thisBf0(i)= bfFun(nullData{:}); % Collect BF under H0
            end
            h1BfAll(j,:) = thisBf1;
            h0BfAll(j,:) = thisBf0;
        end
        results.H1.bf.all= h1BfAll; % Assign outside parfor
        results.H0.bf.all = h0BfAll;
        % Calculate the minimum number of samples needed to reach
        % H1 evidence boundary in p.Restuls.fractionSuccess of the
        % experiments
        ix  = mean(results.H1.bf.all>upperBF)>p.Results.pSuccess; %
        if any(ix)
            results.H1.N.min = Inf;
        else
            results.H1.N.min = min(p.Results.N(ix));
        end
        ix  = mean(results.H0.bf.all<lowerBF)>p.Results.pSuccess; %
        if any(ix)
            results.H0.N.min = Inf;
        else
            results.H0.N.min = min(p.Results.N(ix));
        end
        results.H0.pFalse =  nanmean(results.H0.bf.all>upperBF);  % Upper evidence boundary under H0 = false positives
        results.H1.pFalse =  nanmean(results.H1.bf.all<lowerBF);  % Lower evidence boundary under H1 - false negatives
        results.H0.bf.median  = median(results.H0.bf.all);
        results.H1.bf.median = median(results.H1.bf.all);
        results.H0.pTrue =  nanmean(results.H0.bf.all<lowerBF);  % Lower evidence boundary under H0 = true positives
        results.H1.pTrue =  nanmean(results.H1.bf.all>upperBF);  % Upper evidence boundary under H1 - true positives
        
    end
    results.parms = p.Results; % Store input parser results
end


% If requested show graphs
if showGraph
    
    ticks= ([1/10 1/3 1  3 10 30 10.^(2:8)] );
    tickLabels = {'1/10','1/3', '1',' 3','10','30','100','1000','10^4','10^5','10^6','10^7','10^8'};
    clf;
    if p.Results.sequential
        % Show the "trajectories"  - BF as a funciton of sample
        % size , terminating at the boundaries.
        R =2;C=1;
        out = ticks<lowerBF | ticks>upperBF;
        ticks(out)=[];
        tickLabels(out) =[];
        for i=0:1
            subplot(R,C,1+i);
            if i==0
                toPlot = results.H1.bf.all;
                hyp ='H1';
                if isa(p.Results.effectSize,'function_handle')
                    es = func2str(p.Results.effectSize);
                else
                    es = num2str(p.Results.effectSize);
                end
            else
                toPlot = results.H0.bf.all;
                hyp = 'H0';
                es = '0';
            end
            toPlot(toPlot>upperBF) = upperBF;
            toPlot(toPlot<lowerBF) = lowerBF;
            if size(toPlot,1)>100
                % Limit to 100 trajectories
                toPlot100 = toPlot(randsample(size(toPlot,1),100),:);
            else
                toPlot100 = toPlot;
            end
            plot(p.Results.N,toPlot100','Color',0.5*ones(1,3))
            hold on
            plot([1 maxN],upperBF*ones(1,2),'k--')
            plot([1 maxN],lowerBF*ones(1,2),'k--')
            plot(repmat([1 maxN]',[1 numel(ticks)]), repmat(ticks,[2 1]),'k:')
            set(gca,'YScale','Log','YLim',[lowerBF*0.7 upperBF*1.3],'YTick',ticks,'YTickLabels',tickLabels)
            xlabel 'Sample Size'
            ylabel 'Bayes Factor (BF_{10})'
            title (['BFDA (' test '): Under ' hyp ': Effect Size= ' es ' (nrMC = ' num2str(p.Results.nrMC) ')']);
            h1Boundary  = mean(max(toPlot,[],2)>=upperBF);
            text(min(p.Results.N),1.1*upperBF,sprintf('%d%% stopped at H1 Boundary',round(100*h1Boundary)),'FontWeight','Bold');
            h0Boundary = mean(min(toPlot,[],2)<=lowerBF);
            text(min(p.Results.N),0.9*lowerBF,sprintf('%d%% stopped at H0 Boundary',round(100*h0Boundary)),'FontWeight','Bold');
        end
        
    else
        % 1 or more fixed N designs.
        % Show evidence histograms for each N (columns) under H1 (top
        % row) and under H0 (bottom row).
        C =nrN;
        R =2;
        nrBins = max(10,p.Results.nrMC/50);
        maxlBf = round(log10(prctile(results.H1.bf.all(:),99)));
        minlBf = round(log10(prctile(results.H0.bf.all(:),1)));
        bins= linspace(minlBf,maxlBf,nrBins);
        out = log10(ticks)<minlBf | log10(ticks)>maxlBf;
        
        ticks(out)=[];
        tickLabels(out) =[];
        for j=0:1
            if j==0
                toPlot = results.H1.bf.all;
                hyp ='H1';
                if isa(p.Results.effectSize,'function_handle')
                    es = func2str(p.Results.effectSize);
                else
                    es = num2str(p.Results.effectSize);
                end
            else
                toPlot = results.H0.bf.all;
                hyp = 'H0';
                es = '0';
            end
            for i=1:nrN
                subplot(R,C,C*j+i);
                [n,x] = hist(log10(toPlot(:,i)),bins);
                plot(10.^x,n./sum(n))
                set(gca,'XTick',ticks,'XTickLabel',tickLabels,'XScale','Log')
                hold on
                xlabel 'Bayes Factor (BF_{10})'
                ylabel 'Density'
                title([ hyp ': N = ' num2str(p.Results.N(i)) ', Effect Size = ' es]);
            end
        end
        % Equal x-axes
        allAx = get(gcf,'Children');
        xlims = get(allAx,'XLim');
        xlims = cat(1,xlims{:});
        set(allAx,'XLim',[round(min(xlims(:,1))) round(max(xlims(:,2)))]);
        %         h = annotation(gcf,'textbox',[0.3 0.95 0.4 0.025],'String',...
        %             ['BFDA: ' test ' (nrMC = ' num2str(p.Results.nrMC) ') Minimal N for ' num2str(p.Results.pSuccess*100) '% success = ' num2str(results.H1.N.min) ]);
        h.HorizontalAlignment = 'Center';
        for i=allAx(:)'
            % Show evidence boundary
            plot(i,repmat(p.Results.evidenceBoundary,[2 1]),repmat(ylim(i)',[1 size(p.Results.evidenceBoundary,2)]),'k--')
        end
    end
end
end
