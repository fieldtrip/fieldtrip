%% specification and inference for an unrolled discrete dbn

clear all

rand('twister',1); randn('state',1);

n = 5; % number of nodes per slice
nslices = 5; 
nobservations = 5; % number of random observations

query = 4 + (nslices-1)*n; % query node 3 in last slice

% create some evidence
evidence = nan(nslices,n);
x = randperm(numel(evidence));
for j=1:nobservations
    evidence(x(j)) = round(rand) + 1;
end

dbn = dbnet.random(n,'continuous',[]);

% test complete learning;
% data is a sequence of observations for all nodes
data = rand(100,n);
data(data > 0.5) = 2;
data(data <= 0.5) = 1;

% learn parameters which are coupled between slices
dbn = dbn.learn_parameters(data);

% now the dbn becomes a standard bayesian network
dbn = dbn.unroll(nslices);

% standard inference engine
ie = discrete_jtree_ie(dbn);

% add evidence to the inference engine
ie.enter_evidence(evidence);

% compute a marginal
m = normalize(ie.marginalize(query))

