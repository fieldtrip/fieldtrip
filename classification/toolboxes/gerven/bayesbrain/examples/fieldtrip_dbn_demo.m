%% Using dynamic Bayesian networks together with FieldTrip data
% This example demonstrates how to use neuroimaging data obtained from FieldTrip
% together with the Bayesbrain toolbox. In the example, we make use of covert
% attention data of one subject that has already been frequency analyzed.
% Note that this is trial based data so we can build finite size models
% that capture behaviour within a trial.
%
% The data consists of 7 different frequencies at 274 channels at time
% points [-0.5 0 0.5 1 1.5 2 2.5]. We can expect evoked response after the
% cue and alpha modulation after about 1 second
%
% In this particular example, we will construct a dynamic Bayesian network 
% that captures the evolution over time of the random variables.
%
% Copyright (C) 2008  Marcel van Gerven
%

%% Create and use a DBN
function fieldtrip_dbn_demo()

%%
% Load frequency data and convert the data to a format that can be used by Bayesbrain.
% A Bayesian network is a static model, so we will take out the time data
% and focus only on the 12 Hz band in two channels in left and right hemisphere.

load covattfrq1

% left and right channels
l = find(ismember(left.label,'MLO32'));
r = find(ismember(left.label,'MRO32'));

% We take the log, center, and scale the data to make it better behaved
% Note that we add an arbitrary variable which will accommodate the discrete
% variable
left = squeeze(log(left.powspctrm(:,[l r 1],3,:)));
right = squeeze(log(right.powspctrm(:,[l r 1],3,:)));

% Reshape so each time slice is concatenated after one another
sz = size(left);
left = reshape(left,[sz(1) prod(sz(2:end))]);
sz = size(right);
right = reshape(right,[sz(1) prod(sz(2:end))]);

% Each third variable is the discrete node 
left(:,3:3:size(left,2)) = 1;
right(:,3:3:size(right,2)) = 2;

%%
% Now we can create a very simple model that consists of one discrete
% parent (the attention condition) and two continuous children per slice
data = [left; right];
clear left; clear right;

%%
%  Create a two-slice DBN where nodes in the second slice have
% auto-connections
factors = cell(1,6);
factors{1} = gaussian_cpd(1,[],3,[0; 0],{[]; []},[1; 1]);
factors{2} = gaussian_cpd(2,[],3,[0; 0],{[]; []},[1; 1]);
factors{3} = multinomial_cpd(3,[],[0.5; 0.5]);
factors{4} = gaussian_cpd(4,1,6,[0; 0],{0; 0},[1; 1]);
factors{5} = gaussian_cpd(5,2,6,[0; 0],{0; 0},[1; 1]);
factors{6} = multinomial_cpd(6,3,[0.5 0.5; 0.5 0.5]);

% optionally add names to the factors
factors{1}.name = 'MLO32';
factors{2}.name = 'MRO32';
factors{3}.name = 'orientation';
factors{3}.statenames = {'left attention' 'right attention'};
factors{4}.name = 'MLO32';
factors{5}.name = 'MRO32';
factors{6}.name = 'orientation';
factors{6}.statenames = {'left attention' 'right attention'};

%%
% Create simple dynamic bayes net
dbn = dbnet(factors);

%% 
% Unroll it to accommodate all 7 slices
dbn = dbn.unroll(7,-0.5:0.5:2.5);

%%
% Show the unrolled DBN

dbn.write('tmpdbn','dot','extension','ps');

%% 
% This is what the plot would look like
%
% <<tmpdbn.jpg>>

%%
% Treat DBN as BN
bn = bayesnet(dbn);

%%
% Learn parameters which are coupled between slices
bn = bn.learn_parameters(data);

%% 
% plot parameter uncertainty for the first random variable

mu = bn.factors{1}.ess.mu.value{1};
sigma = sqrt(bn.factors{1}.sigma2(1)/bn.factors{1}.ess.tau.value{1});
subplot(1,2,1);
fplot(@(x)(normpdf(x,mu,sigma)),[mu-3*sigma mu+3*sigma]);
title('uncertainty in the mean');

a = bn.factors{1}.ess.rho.value{1}/2;
b = bn.factors{1}.ess.phi.value{1}/2;
s = bn.factors{1}.sigma2(1);
subplot(1,2,2);
fplot(@(x)((b^a/gamma(a)) * x^(-a-1)* exp(-b/x)),[s - s/2 s + s/2])
title('uncertainty in the variance');

%% 
% show likelihood

bn.loglik(data)

%% 
% compare without assuming that parameters are coupled between slices
dbn = dbnet(factors,'coupled',false);
bn2 = bayesnet(dbn.unroll(7,-0.5:0.5:2.5));

%%
% show evolution of the mean and variance estimates for one variable
bn2 = bn2.learn_parameters(data);

%% 
% likelihood is a bit better without the coupling
bn2.loglik(data)

%% 
% let's check the evolution of the means and standard deviations

mus = zeros(2,2,7);
sds = zeros(2,2,7);
for s=1:2
    for j=1:7   
        mus(s,1,j) = bn2.factors{(j-1)*3+s}.mu(1);
        mus(s,2,j) = bn2.factors{(j-1)*3+s}.mu(2);        
        sds(s,1,j) = sqrt(bn2.factors{(j-1)*3+s}.sigma2(1));
        sds(s,2,j) = sqrt(bn2.factors{(j-1)*3+s}.sigma2(2));        
    end
end

figure
subplot(1,2,1);
errorbar(-0.5:0.5:2.5,squeeze(mus(1,1,:)),squeeze(sds(1,1,:)));
hold on;
errorbar(-0.5:0.5:2.5,squeeze(mus(1,2,:)),squeeze(sds(1,2,:)),'r');
title('MLO32');
legend('left attention','right attention');

subplot(1,2,2);
errorbar(-0.5:0.5:2.5,squeeze(mus(2,1,:)),squeeze(sds(2,1,:)));
hold on;
errorbar(-0.5:0.5:2.5,squeeze(mus(2,2,:)),squeeze(sds(2,2,:)),'r');
title('MRO32');
legend('left attention','right attention');

%% 
% let's see how the variables depend on their past state

betas = zeros(2,2,6);
for s=1:2
    for j=1:6   
        betas(s,1,j) = bn2.factors{j*3+s}.beta{1};
        betas(s,2,j) = bn2.factors{j*3+s}.beta{2};        
    end
end

figure
subplot(1,2,1);
plot(0:0.5:2.5,squeeze(betas(1,1,:)));
hold on;
plot(0:0.5:2.5,squeeze(betas(1,2,:)),'r');
title('MLO32');
legend('left attention','right attention');

subplot(1,2,2);
plot(0:0.5:2.5,squeeze(betas(2,1,:)));
hold on;
plot(0:0.5:2.5,squeeze(betas(2,2,:)),'r');
title('MRO32');
legend('left attention','right attention');

