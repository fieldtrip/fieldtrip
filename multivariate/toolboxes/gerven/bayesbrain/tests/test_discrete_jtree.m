n = 150;

rand('twister',1); randn('state',1);

bn = bayesnet.random(n,'continuous',[]);

ie = discrete_jtree_ie(bn);

evid = nan(1,n); evid(1) = 1; evid(7) = 2;

ie.enter_evidence(evid);

m = normalize(ie.marginalize(4))