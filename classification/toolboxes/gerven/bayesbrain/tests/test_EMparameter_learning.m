% demonstrates learning

rand('twister',1); randn('state',1);

factors = cell(1,1);
factors{1} = gaussian_cpd(1,[],[],-10,{[]},1);

bn = bayesnet(factors);

% simulate data
data = bn.simulate(300);
data = [nan(size(data,1),1) data];

% now try to reconstruct with a different model
factors = cell(1,2);
factors{1} = multinomial_cpd(1,[],[0.5 0.5]');
factors{2} = gaussian_cpd(2,[],1,[1; 1],cell(2,1),[1; 1]);
bn = bayesnet(factors);

bn.emiter = 2;
bn = bn.learn_parameters(data);

bn.factors{2}.mu
bn.factors{2}.sigma2

% % demonstrates learning
% 
% factors = cell(1,2);
% factors{1} = multinomial_cpd(1,[],[0.5 0.5]');
% factors{2} = gaussian_cpd(2,[],1,[-10; 10],cell(2,1),[1; 1]);
% 
% bn = bayesnet(factors);
% 
% % simulate data
% data = bn.simulate(1000);
% 
% % now try to reconstruct
% 
% data(:,1) = nan;
% 
% bn.emiter = 10;
% bn = bn.learn_parameters(data);
% 
% bn.factors{2}.mu
% bn.factors{2}.sigma2

% % demonstrates learning
% 
% factors = cell(1,3);
% factors{1} = multinomial_cpd(1,[],[0.5 0.5]');
% factors{2} = multinomial_cpd(2,[],[0.5 0.5]');
% factors{3} = gaussian_cpd(3,[],[1 2],[-10 10; -20 20],cell(2,2),[1 1; 1 1]);
% 
% bn = bayesnet(factors);
% 
% % simulate data
% data = bn.simulate(1000);
% 
% % now try to reconstruct
% 
% data(:,1) = nan;
% 
% bn.emiter = 10;
% bn = bn.learn_parameters(data);
% 
% bn.factors{3}.mu
% bn.factors{3}.sigma2

% n = 10; %8
% 
% rand('twister',1); randn('state',1);
% 
% bn = bayesnet.random(n,'continuous',[2 4 6 8]);
% 
% % test complete learning
% 
% data = rand(100,10);
% data(data > 0.5) = 2;
% data(data <= 0.5) = 1;
% 
% % create unobserved data
% x = randperm(numel(data));
% data(x(1:10)) = nan;
% 
% bn = bn.learn_parameters(data);
% 
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