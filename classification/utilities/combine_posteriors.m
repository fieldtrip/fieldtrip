function post = combine_posteriors(cpost,method)
%COMBINE_POSTERIORS combines multiple posteriors using some combination
% method.
%
%   post = combine_posteriors(cpost,method)
%
%   cpost is a cell array of posteriors
%   method is one of
%   'product' : normalized product of posteriors (default)
%   'majority' : majority vote
%   'mean' : normalized sum of posteriors
%   'concatenate' : concatenate posteriors to be used as data
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: combine_posteriors.m,v $
%

% return if cpost is just one posterior
if iscell(cpost) && length(cpost) == 1
    post = cpost{1};
    return
elseif ~iscell(cpost)
    post = cpost;
    return
end

if nargin < 2
    method = 'product';
end

% now combine the posteriors using some combination rule
switch lower(method)

    case 'majority' % majority vote
        
        post = zeros(size(cpost{1}));
        for k=1:length(cpost)
            [temp,pcls] = max(cpost{k},[],2);
            for p=1:length(pcls)
                post(p,pcls(p)) = post(p,pcls(p)) + 1;
            end
        end

        % resolve ties
        for p=1:size(post,1)
            m = find(ismember(post(p,:),max(post(p,:))));
            post(p,:) = 0;
            if length(m) > 1
                m = m(ceil(rand*length(m)));
            end
            post(p,m) = 1;
        end
        
        return % no normalization required
       
    case 'product' % normalized product of probabilities
        
        post = ones(size(cpost{1}));
        for k=1:length(cpost)
            post = post .* cpost{k};
        end
        
    case 'mean' % normalized sum of probabilities
        
        post = zeros(size(cpost{1}));
        for k=1:length(cpost)
            post = post + cpost{k};
        end
        
    case 'concatenate' % concatenate nclasses-1 probabilities to act as input to another classifier

        nclasses = size(cpost{1},2) - 1;

        post = zeros(size(cpost{1},1),nclasses*length(cpost));
        for k=1:length(cpost)
            post(:,(k-1)*nclasses + (1:nclasses)) = cpost{k}(:,1:nclasses);
        end

        return % no normalization required
        
    otherwise 
        error('unknown option for COMBINE_POSTERIORS');

end

 % normalize
 post = post ./repmat(sum(post,2),[1 size(post,2)]);
       

