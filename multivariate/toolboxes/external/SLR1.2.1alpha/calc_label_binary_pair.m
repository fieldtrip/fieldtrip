function	[label_est, Prob] = calc_label_binary_pair(X, w, pairs, mode)
% Calculate label for inputs from one-versus-one (pair-wise) classifiers  
%
% [label_est, Prob] = calc_label_binary_pair(X, w, pair, mode)
%
% --- Input 
% X(:,t) : input data [Xdim x Ndata]
% w(:,n) : n-th binary model
% pairs(n,:) : n-th pair 
%
% --- Optional Input
% mode [0] : mode for combination 
%      = 0 : multiplicative combination of M binary models [default]
%      = 1 : summation combination of M binary models 
%
% --- Output
% label_est(t) : class ID for X(:,t) [1 x Ndata] (id = {1:M})
% Prob(m,t) : m-th class prbability for X(:,t) [M x Ndata]
%         = P(Zall(m)=1|X(t))
%
% 2009/06/05 OY
% * revised for SLR toolbox
% 2008/12/30 Masa-aki Sato
% * original version
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


if nargin < 4, mode = 0; end;

MinProb = eps; % = 0.0001;

Nclass = max(pairs(:));   % # of classes
Nmodel = size(pairs,1);  % # of pairs

[Ndata, Ndim] = size(X);

% Binary probability for each model

switch mode
case	0
	logP = zeros(Ndata,Nclass);
case	1
	P = zeros(Ndata,Nclass);
end

% combine Nmodel binary models
for nn=1:Nmodel
    % set to prior probability (assume balanced samples)
    px = 1/Nclass*ones(Ndata,Nclass); 
    
	% probability of binary models
	[tmp_id, pt] = calc_label(X, w(:,nn));  % pt : probability observing label c2 [1*Ndata] 
	
    pair = pairs(nn,:);
    
    pxt = pt(:,1);
	pxt = min( max(pxt, MinProb), 1-MinProb);% 0 < pxt < 1
	
	% If pxt = 1/2, px(m) = 1/Nclass : No information
	px(:,pair(1)) = (1-pxt)*(2/Nclass);
	px(:,pair(2)) = pxt*(2/Nclass);
	
	switch mode
	case	0
		logP = log(px) + logP;
	case	1
		P = P + px;
	end
end

switch mode
    case	0 
        Prob = smlog2prob(logP);                  
    case	1   % normalization
        Prob = P ./ repmat(sum(P,2), [1,Nclass]);
end

[pmax, label_est ] = max(Prob,[],2);



