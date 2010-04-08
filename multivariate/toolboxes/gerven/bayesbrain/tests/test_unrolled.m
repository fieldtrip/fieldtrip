% demonstrates inference for a random dynamic Bayesian network

n = 5; %8
nslices = 5;
nobs = 5;

rand('twister',3); randn('state',1);

query = ceil(n*rand); % query node

% create evidence
evidence = nan(nslices,n);
x = randperm(numel(evidence));
for j=1:nobs
    evidence(x(j)) = round(rand) + 1;
end

dbn = dbnet.random(n,'continuous',[]);

% test complete learning
data = rand(100,n);
data(data > 0.5) = 2;
data(data <= 0.5) = 1;

dbn = dbn.learn_parameters(data);

if 0
   query = query + (nslices-1)*n; % query node
   dbn = dbn.unroll(nslices);
   ie = discrete_jtree_ie(dbn);
else
    ie = filtering_ie(dbn);
end

% add evidence to the inference engine
ie.enter_evidence(evidence);

% compute a marginal
m = normalize(ie.marginalize(query))

