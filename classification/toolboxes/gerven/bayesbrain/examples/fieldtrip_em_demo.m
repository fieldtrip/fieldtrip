%% Using the EM algorithm together with FieldTrip data
% This example demonstrates how to use neuroimaging data obtained from FieldTrip
% together with the neurogm toolbox. In the example, we make use of covert
% attention data of one subject that has already been frequency analyzed.
% Note that this is trial based data so we can build finite size models
% that capture behaviour within a trial.
%
% The data consists of 7 different frequencies at 274 channels at time
% points [-0.5 0 0.5 1 1.5 2 2.5]. We can expect evoked response after the
% cue and alpha modulation after about 1 second
%
% In this particular example, we will construct a standard Bayesian network and
% demonstrate the use of the EM algorithm to get better estimates of the
% distributions.
%
% Copyright (C) 2008  Marcel van Gerven
%

%% Compare log likelihood of a Bayesian network before/after learning
function fieldtrip_em_demo()

%%
% Load frequency data and convert the data to a format that can be used by NEUROGM.
% A Bayesian network is a static model, so we will take out the time data
% and focus only on the 12 Hz band in two channels in left and right
% hemisphere.

load freqli; % left attention
load freqri; % right attention

% left and right channels
l = find(ismember(freqLI.label,'MLO32'));
r = find(ismember(freqRI.label,'MRO32'));

% We take the log and ztransform to make the data better behaved. Essential
% for the EM algorithm!
datal = zscore(log((squeeze(nan_mean(freqLI.powspctrm(:,[l r],3,:),4)))));
datar = zscore(log((squeeze(nan_mean(freqRI.powspctrm(:,[l r],3,:),4)))));
clear freqli; clear freqri;

%%
% Now we can create a very simple model that consists of one discrete
% parent (the attention condition) and two continuous children who each
% have an extra discrete variable which allows more complex distributions
% to be estimated.
data = [[datal ones(size(datal,1),1) nan(size(datal,1),2) ]; [datar 2*ones(size(datar,1),1) nan(size(datar,1),2)]];

%%
%  Create the random variables; they should follow the data ordering
factors = cell(1,5);
factors{1} = gaussian_cpd(1,[],[3 4],zeros(2,2),cell(2,2),100*ones(2,2));
factors{2} = gaussian_cpd(2,[],[3 5],zeros(2,2),cell(2,2),100*ones(2,2));
factors{3} = multinomial_cpd(3,[],[0.5; 0.5]);
factors{4} = multinomial_cpd(4,[],rand(2,1));
factors{5} = multinomial_cpd(5,[],rand(2,1));

% optionally add names to the factors
factors{1}.name = 'MLO32';
factors{2}.name = 'MRO32';
factors{3}.name = 'orientation';
factors{3}.statenames = {'left attention' 'right attention'};
factors{4}.name = 'latent1';
factors{5}.name = 'latent2';

%%
% Create simple bayes net
bn = bayesnet(factors);

%%
% Write graph structure to .ps file (requires installation of GraphViz
% library)
bn.write('tmpem','dot','extension','ps');

%% 
% This is what the plot would look like
%
% <<tmpem.jpg>>

%% 
% Learn parameters from incomplete data (invoke EM algorithm)

bn = bn.learn_parameters(data);

%%
% Plot of the  mixture distributions

% variable of interest
obj = bn.factors{2};

% mixture components and histograms
mix = bn.factors{5};

lim = [min(obj.mu(:)) - 3*sqrt(max(abs(obj.sigma2(:)))) max(obj.mu) + 3*sqrt(max(abs(obj.sigma2(:))))];
subplot(1,2,1);
h1 = histc(datal(:,2),-2:0.25:2);
bar(-2:0.25:2,h1./(5*max(h1)),'histc');
hold on;
fplot(@(x)(mix.p(1)*normpdf(x,obj.mu(1),sqrt(obj.sigma2(1))) + mix.p(2)*normpdf(x,obj.mu(3),sqrt(obj.sigma2(3)))),lim,'k');
title('MRO32 left attention');

subplot(1,2,2);
h2 = histc(datar(:,2),-2:0.25:2);
bar(-2:0.25:2,h2./(5*max(h2)),'histc');
hold on;
fplot(@(x)(mix.p(1)*normpdf(x,obj.mu(2),sqrt(obj.sigma2(2))) + mix.p(2)*normpdf(x,obj.mu(4),sqrt(obj.sigma2(4)))),lim,'k');
title('MRO32 right attention');

end
