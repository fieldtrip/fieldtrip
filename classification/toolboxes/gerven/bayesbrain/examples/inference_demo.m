%% inference for a discrete Bayesian network

% create random Bayesian network with 100 discrete nodes
n = 100;
bn = bayesnet.random(n,'continuous',[]);

% create junction tree algorithm for discrete networks
% we could also use hugin_ie or exhaustive_ie
ie = discrete_jtree_ie(bn);

% construct evidence (nan = unobserved)
evid = nan(1,n);
evid(1) = 1; evid(7) = 2;

% enter evidence into the inference engine
ie.enter_evidence(evid);

% compute marginal for node 4
m = normalize(ie.marginalize(4))

%% inference for a hybrid Bayesian network

n = 8;

rand('twister',1); randn('state',1);

bn = bayesnet.random(n,'continuous',[1 3 5 6]);

% we could also use hugin_ie or exhaustive_ie
ie = canonical_jtree_ie(bn);

% add evidence
evid = nan(1,n); evid(1) = 1; evid(5) = 2;

ie.enter_evidence(evid);

m = normalize(ie.marginalize(3));

