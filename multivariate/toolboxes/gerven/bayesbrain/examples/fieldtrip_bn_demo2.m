%% Using Bayesian networks together with FieldTrip data
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
% demonstrate its use. However, we will represent the power as 
% power / mean power per channel-frequency pair
%
% Copyright (C) 2008  Marcel van Gerven
%

%% try alternative representation of power
function fieldtrip_bn_demo2()

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

% We take the log to make the data better behaved
datal = (squeeze(nan_mean(freqLI.powspctrm(:,[l r],3,:),4)));
datar = (squeeze(nan_mean(freqRI.powspctrm(:,[l r],3,:),4)));
clear freqli; clear freqri;

%%
% Now we can create a very simple model that consists of one discrete
% parent (the attention condition) and two continuous children
data = [datal; datar];
data = data./repmat(mean(data,1),[size(data,1) 1]);
data = [data [ones(size(datal,1),1); 2*ones(size(datar,1),1)]];
clear datal; clear datar;

%%
%  Create the random variables; they should follow the data ordering
factors = cell(1,3);
factors{1} = gaussian_cpd(1,[],3,[0; 0],{[]; []},[1; 1]);
factors{2} = gaussian_cpd(2,[],3,[0; 0],{[]; []},[1; 1]);
factors{3} = multinomial_cpd(3,[],[0.5; 0.5]);

% optionally add names to the factors
factors{1}.name = 'MLO32';
factors{2}.name = 'MRO32';
factors{3}.name = 'orientation';
factors{3}.statenames = {'left attention' 'right attention'};

%%
% Create simple bayes net
bn = bayesnet(factors);

%%
% Write graph structure to .ps file (requires installation of GraphViz
% library)
bn.write('tmpbn','dot','extension','ps');

%% 
% This is what the plot would look like
%
% <<tmp.jpg>>

%%
% Learn parameters from complete data
bn = bn.learn_parameters(data);

%%
% Log likelihood
bn.loglik(data)

%%
% Plot the estimated prior distributions with continuous ones of the form
%
% $$ \mathcal{N}(x; \mu_c, \sigma_c)$$
%
        
subplot(1,3,1);
bn.factors{1}.plot();
legend('left attention','right attention');
subplot(1,3,2);
bn.factors{2}.plot();
legend('left attention','right attention');
subplot(1,3,3);
bn.factors{3}.plot();
set(gcf,'Position',[100 100 1500 400]);

%%
% Create an inference engine

ie = canonical_jtree_ie(bn);

%% 
% Add some evidence

ie.enter_evidence([nan 1.5 nan]);

%%
% Compute marginals

m1 = normalize(ie.marginalize(1));
m3 = normalize(ie.marginalize(3));

%% 
% Plot the marginals after evidence propagation

figure
subplot(1,2,1);
m1.plot();
subplot(1,2,2);
m3.plot();
