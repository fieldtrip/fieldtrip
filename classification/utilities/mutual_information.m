function mi = mutual_information(data,design)
% MUTUAL_INFORMATION compute mutual information between features and class, 
% assuming that feature values are normally distributed conditional on the class
%
%   mi = mutual_information(data,design)
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: mutual_information.m,v $
%
    % number of class labels
    csize = max(design(:,1));

    % compute prior for class variable
    pc = zeros(1,csize);
    for k=1:csize
        pc(k) = sum(design(:,1) == k);
    end
    pc = pc ./ size(data,1);

    % precompute indices
    di = cell(1,csize);
    for k=1:csize, di{k} = (design(:,1) == k); end

    mi = zeros(1,size(data,2));

    for m=1:size(data,2)

        % use numerical integration to approximate p(x) log p(x)

        for k=1:csize
            mi(m) = mi(m) - 0.5 * pc(k) * log2(var(data(di{k},m)));
        end

        mi(m) = mi(m) - 0.5 * log2(2*pi*exp(1));

        % construct p(x) log p(x)  = (\sum_c p(c) p(x|c)) log (\sum_c p(c) p(x|c))

        % create \sum_c p(c) p(x|c) anonymous function
        f = @(x,k)( pc(k) .*  mynormpdf(x,nanmean(data(di{k},m)),mynanstd(data(di{k},m))));
        g = @(x)(sum(cell2mat(transpose(arrayfun(@(k)(f(x,k)),1:csize,'UniformOutput',false))),1));
        h = @(x)(pp(g(x)));

        mi(m) = mi(m) - quadgk(h,-Inf,Inf);
    end
end

function p = mynormpdf(x,mu,sigma)
% MYNORMPDF implements the normal distribution with mean mu and standard
% deviation sigma evaluated at x
%
% p = mynormpdf(x,mu,sigma)
%

    if sigma
        p = 1./(sqrt(2.*pi).*sigma) .* exp(- (x-mu).^2/(2*sigma.^2));
    else
        p = zeros(1,length(x));
        p(x==mu) = 1;
    end
end

function y = pp(x)

    x(~x(:)) = 1;
    y = x .* log2(x);
    y(isnan(y(:))) = 0;
end
