function mi = conditional_mutual_information(data,design,conditional,lim)
% CONDITIONAL_MUTUAL_INFORMATION compute mutual information between features 
% and class conditioned on another feature, assuming that feature values are 
% normally distributed conditional on the class
%
%   mi = conditional_mutual_information(data,design,conditional,lim)
%
%   data is data for features without conditional
%   design is design matrix containing class labels
%   conditional is a column of data for the conditional
%   lim is the limit of the rectangle at which to evaluate the integral
%
%   TODO:
%   This code is not yet working!
%
%   EXAMPLE:
%   
%   dat = rand(20,1); conditional = rand(20,1); design = [ones(size(dat,1)/2,1); 2*ones(size(dat,1)/2,1)];
%   conditional_mutual_information(dat,design,conditional)
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: conditional_mutual_information.m,v $
%

    % limits of the rectangle in which to evaluate integral
    % data should be centered around zero
    if nargin < 4
        lim = 1;
    end

    % number of class labels
    csize = max(design(:,1));

    % compute prior for class variable
    pc = zeros(1,csize);
    for k=1:csize, pc(k) = sum(design(:,1) == k); end
    pc = pc ./ size(data,1);

    % precompute indices
    di = cell(1,csize);
    for k=1:csize, di{k} = (design(:,1) == k); end
           
    mi = zeros(1,size(data,2));
    
    for m=1:size(data,2)

        % create mu_c and Sigma_c
        mu = cell(1,csize);
        Sigma = cell(1,csize);
        for k=1:csize
            mu{k} = [nanmean(data(di{k},m)); nanmean(conditional(di{k}))];
            Sigma{k} = nancov([data(di{k},m) conditional(di{k})]);
        end
        
        % use numerical integration to approximate p(x) log p(x)

        for k=1:csize
            mi(m) = mi(m) - pc(k) * log2(sqrt((2*pi*exp(1))^2*det(Sigma{k})));            
        end
        
        % numerical approximation of second term
                
        for k=1:csize
        
            % construct function
            f2 = @(x,y)(pp(mynormpdf([x; y*ones(1,size(x,2))],mu{k},Sigma{k}),mynormpdf(y,mu{k}(2),Sigma{k}(4))));

            tlim = lim*2;
            intgrl = nan;
            while isnan(intgrl) || isinf(intgrl)
                tlim = tlim/2;
                intgrl = dblquad(f2,-tlim,tlim,-tlim,tlim);
            end
            
            mi(m) = mi(m) - pc(k) * intgrl;
        end
        
        % numerical approximation of third term        
        for k=1:csize
           
            g1 = @(x,y,k)( pc(k) .*  mynormpdf([x; y*ones(1,size(x,2))],mu{k} + (Sigma{k}(3)/Sigma{k}(4))...
                *(y - mu{k}(2)),Sigma{k}(1) - (Sigma{k}(3).^2/Sigma{k}(4))) );
            g2 = @(x,y)(sum(cell2mat(transpose(arrayfun(@(k)(g1(x,y,k)),1:csize,'UniformOutput',false))),1));
            g3 = @(x,y)(pp(mynormpdf([x; y*ones(1,size(x,2))],mu{k},Sigma{k}),g2(x,y)));
            
            tlim = lim*2;
            intgrl = nan;
            while isnan(intgrl) || isinf(intgrl)
                tlim = tlim/2;
                intgrl = dblquad(g3,-tlim,tlim,-tlim,tlim);
            end
            
            mi(m) = mi(m) - pc(k) * intgrl;
        end

    end
end

function p = mynormpdf(x,mu,sigma)
% MYNORMPDF implements the normal distribution with mean mu and standard
% deviation sigma evaluated at x
%
% p = mynormpdf(x,mu,sigma)
%

    if isscalar(x)
        p = 1./(sqrt(2.*pi).*sigma) .* exp(- (x-mu).^2/(2*sigma.^2));
    else

        n = size(x,1);
        
        p = zeros(1,size(x,2));

        for j=1:size(x,2)

            p(j) = (1./(((2*pi)^(n/2))*sqrt(det(sigma)))) * exp(-0.5*(x(:,j)-mu)'*inv(sigma)*(x(:,j)-mu));
        end

    end
    
    p(isnan(p)) = 0;
    
end

function y = pp(a,b)

    y = a .* log2(b);    
    y(isnan(y)) = 0;
end