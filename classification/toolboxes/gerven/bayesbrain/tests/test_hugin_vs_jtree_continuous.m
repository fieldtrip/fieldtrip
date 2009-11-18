% demonstrates inference for a random Bayesian network

% n = 20; %8
% 
% rand('twister',x); randn('state',1);
% 
% bn = bayesnet.random(n,'continuous',[2 3 4]);
% 
% ie = hugin_ie(bn);
% 
% evid = nan(1,n); evid(2) = 1; evid(3) = 2;
% 
% %for i=1:50
% tic
%     ie.enter_evidence(evid);
% %end
% 
% m = normalize(ie.marginalize(4));
% toc
% ie2 = hugin_ie(bn);
% 
% ie2.enter_evidence(evid);
% 
% m2 = normalize(ie2.marginalize(4));
% 
% if ~all(abs(m.mu - m2.mu) < 0.01)
%     error('different');
% end


n = 20; %8

rand('twister',x); randn('state',1);

bn = bayesnet.random(n,'continuous',[1 3 5 6]);

ie = canonical_jtree_ie(bn);

evid = nan(1,n); evid(1) = 1; evid(5) = 2;

%profile on
%tic
%for i=1:50
ie.enter_evidence(evid);
%end
%toc
%profile off
%profile report

m = normalize(ie.marginalize(3));

% ie2 = hugin_ie(bn);
% 
% ie2.enter_evidence(evid);
% 
% m2 = normalize(ie2.marginalize(3));
% 
% 
% if ~all(abs(m.mu - m2.mu) < 0.01)
%     error('different');
% end
% 
% 
% n = 20; %8
% 
% rand('twister',x); randn('state',1);
% 
% bn = bayesnet.random(n,'continuous',[1 3 5 6]);
% 
% ie = canonical_jtree_ie(bn);
% 
% evid = nan(1,n); evid(1) = 1; evid(5) = 2;
% 
% ie.enter_evidence(evid);
% 
% m = normalize(ie.marginalize(3));
% 
% ie2 = hugin_ie(bn);
% 
% ie2.enter_evidence(evid);
% 
% m2 = normalize(ie2.marginalize(3));
% 
% 
% if ~all(abs(m.mu - m2.mu) < 0.01)
%     error('different');
% end