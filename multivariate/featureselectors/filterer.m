classdef filterer < featureselector
%FILTERER filtering approach to feature selection
%
%   The filterer produces an ordering of the features based on a univariate measure
%   and uses that ordering to greedily determine the 'optimal' feature set
%   or to select the best m features.
%
%   PARAMETERS:
%    criterion; % keeps track of the evaluation metric
%    value; % filter function outputs
%    order; % ordering of the variables
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    nfeatures = Inf; % maximum number of used features
        
    % filter is used to order the features
    %
    % options:
    % [10 100 2 ...] : uses the ordering specified
    % meandiff : order according to the norm of differences of means
    % mutual_information : order according to the mutual information (default)
    % anova : order according to anova test  
    % regression : linearly regress input on output and compute residual
    % correlation : correlation between input and output
    filter = 'anova';    
    
  end
  
  methods
    function obj = filterer(varargin)
      
      obj = obj@featureselector(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
      
      maxfeatures = min(obj.nfeatures,size(X,2));
      
      if obj.verbose
        if isempty(obj.validator)
          fprintf('selecting %d features based on %s filter\n',maxfeatures,obj.filter);
        else
          fprintf('selecting optimal features based on %s filter\n',obj.filter);
        end
      end
        
      % compute function values
      nf = size(X,2);
      p.value = zeros(1,nf);
      for j=1:nf
        
        if obj.verbose
          fprintf('computing %s filter for feature %d of %d\n',obj.filter,j,nf);
        end
        
        switch(obj.filter)
          
          case 'meandiff'
            p.value(j) = filterer.meandiff(X(:,j),Y);
            
          case 'mutual_information'
            p.value(j) = filterer.mutual_information(X(:,j),Y);
            
          case 'anova'
            p.value(j) = filterer.anova(X(:,j),Y);
            
          case 'regression'
            p.value(j) = filterer.regress(X(:,j),Y);
            
          case 'correlation'
            p.value(j) = filterer.correlation(X(:,j),Y);
            
          otherwise
            error('unrecognized filter');
        end
            
      end
      
      % get ordering of the features
      [tmp,p.order] = sort(p.value,'descend');
      
      if isempty(obj.validator)
        % select the best n features
        p.subset = p.order(1:maxfeatures);
      else
        % use procedure to determine the best feature subset
        [p.subset,p.criterion] = obj.select_features(X,Y,p.order);
      end
      
    end
    
    function [subset,criterion] = select_features(obj,X,Y,features)
      % select features
      % assumes availability of evaluation function and feature
      % ordering
      
      subset = [];
      metric = -Inf;
            
      criterion = zeros(1,length(features));
      for f=1:min(obj.nfeatures,length(features))
        
        if obj.verbose
          fprintf('evaluating %d out of %d features; ',f,length(features));
        end
        
        cv = obj.validator.validate(X(:,features(1:f)),Y);
        
        m = cv.evaluate('metric',obj.metric);
        
        if obj.verbose, fprintf('criterion: %.2g\n',m); end
        
        if m > metric
          subset = features(1:f);
          metric = m;
        end
        
        criterion(f) = m;
        
      end
      
    end
    
    function [model,desc] = getmodel(obj)
      % return logical array with ones indicating selected features
      
      indim = obj.indims;
      if length(indims) == 1, indim = [1 indim]; end
      
      m = zeros(indim);
      m(obj.params.subset) = 1;      
      model{1} = m;      
      desc{1} = {'logical array with ones indicating selected features'};
      
      m = zeros(indim);
      m(:) = obj.params.value;
      model{2} = m;
      desc{2} = {sprintf('importance according to filter %s',obj.filter)};
      
    end
    
  end
  
  methods(Static=true)
    
    function dif = meandiff(data,design)
      % MEANDIFF computes the l2 norm of the difference between the means of the data of
      % different classes w.r.t. the global mean per feature
      %
      %   dif = meandiff(data,design)
      
      nclasses = max(design(:,1));

      dif = zeros(nclasses,size(data,2));
      for j=1:nclasses
        dif(j,:) = mynanmean(data(design==j,:));
      end
      dif = dif - repmat(mynanmean(data),[nclasses 1]);
      dif = sqrt(sum(dif.^2));
      
    end
    
    function ip = anova(data,design)
      % ANOVA computes 1-pvalues that data in different classes belong to the
      % same group. I.e., the closer to one, the better.
      %
      %   ip = anova(data,design)
      %
      %   ip = 1 - pvalue
      
      ip = zeros(1,size(data,2));
      for j=1:size(data,2)
        ip(j) = 1 - anova1(data(:,j),design(:,1),'off');
      end
      
    end
    
    function res = regress(X,Y)
      % linearly regress input on output and compute residual

      [a,b,res] = regress(Y,[X ones(size(X,1),1)]);
      
      res = 1./mean(abs(res));
      
    end
    
    function res = correlation(X,Y)
      % compute correlation between input and output

      res = corr(X,Y);
    
    end
    
    function mi = mutual_information(data,design)
      % MUTUAL_INFORMATION compute mutual information between features and class,
      % assuming that feature values are normally distributed conditional on the class
      %
      %   mi = mutual_information(data,design)
      %
      %   Copyright (c) 2008, Marcel van Gerven
      
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
        f = @(x,k)( pc(k) .*  mynormpdf(x,mynanmean(data(di{k},m)),mynanstd(data(di{k},m))));
        g = @(x)(sum(cell2mat(transpose(arrayfun(@(k)(f(x,k)),1:csize,'UniformOutput',false))),1));
        h = @(x)(pp(g(x)));
        
        mi(m) = mi(m) - quadgk(h,-Inf,Inf);
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
      
    end
    
    
    
  end
  
end