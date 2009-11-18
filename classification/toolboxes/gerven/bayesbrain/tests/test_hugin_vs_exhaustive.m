% demonstrates inference for a random Bayesian network

n = 10; %8

rand('twister',x); randn('state',1);

bn = bayesnet.random(n,'continuous',[1 5]);

ie = exhaustive_ie(bn);

evid = nan(1,n); evid(1) = 1; evid(7) = 2;

ie.enter_evidence(evid);

m = normalize(ie.marginalize(4));

ie2 = hugin_ie(bn);

ie2.enter_evidence(evid);

m2 = normalize(ie2.marginalize(4));

if ~all(abs(m.p.value - m2.p.value) < 0.01)
    error('different');
end