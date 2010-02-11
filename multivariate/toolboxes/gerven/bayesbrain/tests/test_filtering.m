% demonstrates inference for a random dynamic Bayesian network

n = 3; %8
nslices = 5;
nobs = 4;

rand('twister',1); randn('state',1);

dbn = dbnet.random(n,'continuous',[]);

% test complete learning
data = rand(100,10);
data(data > 0.5) = 2;
data(data <= 0.5) = 1;

dbn = dbn.learn_parameters(data);

ie = filtering_ie(dbn);

% create evidence
evidence = nan(nslices,n);
x = randperm(numel(evidence));
for j=1:nobs
    evidence(x(j)) = round(rand) + 1;
end

% add evidence to the inference engine
ie.enter_evidence(evidence);

query = ceil(n*rand); % query node

% compute a marginal
m = normalize(ie.marginalize(query))

