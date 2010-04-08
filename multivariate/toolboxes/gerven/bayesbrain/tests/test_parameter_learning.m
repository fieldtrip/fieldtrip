% demonstrates inference for a random dynamic Bayesian network

n = 10; %8

rand('twister',1); randn('state',1);

% HERE!
dbn = dbnet.random(n,'continuous',[]);

% test complete learning

data = rand(100,10);
data(data > 0.5) = 2;
data(data <= 0.5) = 1;

dbn = dbn.learn_parameters(data);

% ie = discrete_jtree_ie(bn);
% 
% evid = nan(1,n); evid(1) = 1; evid(7) = 2;
% 
% ie.enter_evidence(evid);
% 
% m = normalize(ie.marginalize(4));
% 
% ie2 = hugin_ie(bn);
% 
% ie2.enter_evidence(evid);
% 
% m2 = normalize(ie2.marginalize(4));
% 
% if ~all(abs(m.p.value - m2.p.value) < 0.01)
%     error('different');
% end