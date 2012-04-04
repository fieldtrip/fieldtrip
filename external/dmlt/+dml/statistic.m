function s = statistic(stat,D,P)
% STATISTIC test statistics.
%
%   DESCRIPTION
%   Computes test statistic stat (or multiple using a cell array) for a
%   nsamples x nvars design matrix D and a nsamples x nvars prediction
%   matrix P.
%
%   Available statistics:
%
%   'accuracy'    : proportion of correctly classified samples
%   'logprob'     : log probability of the correct class
%   'correlation' : correlation between input and output matrices
%   'contingency' : contingency matrix (rows=true class, columns=predicted class)
%   'confusion'   : confusion matrix (rows=true class, columns=predicted class)
%   'binomial'    : binomial test gives p value for significant
%                   classification; takes uneven class distributions into
%                   account
%   'MAD'         : mean absolute deviation in degrees for angles specified
%                   in radians
%   'RMS'         : root-mean-square error
%   'tlin'        : t-linear statistic to measure goodness of fit for
%                   circular data with angles specified in radians 
%                   using fisher's correlation coefficient.
%   'R2'          : squared correlation; explained variance in the linear regression case
%   'expvar'      : explicit computation of explained variance: 
%                   (var(y) - var(y-yhat)) / var(y)
%                   1 = perfect prediction, <0 = worse than chance
%   'identity'    : computes the proportion of trials that are identified correctly 
%
%   NOTE: notation '-x' with x one of the above is also allowed. This way
%   error measures such as 'MAD' can be used as a performance measure by
%   dml.permutation or dml.gridsearch using '-MAD'
%
%   EXAMPLE
%   s = dml.statistic('accuracy',D,P)
%   s = dml.statistic('-MAD',D,P)
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  if iscell(stat)
    s = cell(length(stat),1);
    for c=1:length(stat)
      s{c} = dml.statistic(stat{c},D,P);
    end
  else
    
    
    minus = stat(1) == '-';
    if minus, stat = stat(2:end); end
    
    nvars = size(D,2);
    nsamples = size(D,1);
    
    switch stat
      
      case 'accuracy'
        
        % compute the total proportion of correctly classified trials
        % note: is insensitive to the relative frequency of
        % examples from each class
        tmp = contingency(D,P);
        s = sum(diag(tmp))./nsamples;
        
      case 'logprob'
        
        % compute summed log probability of the real class
        s = 0;
        for j=1:nsamples
          s = s + log(P(j,D(j)));
        end
        
      case {'correlation','R2'}
        
        s = ones(1,nvars);
        idx = find(~all(D == P));
        for i=1:length(idx)
          s(idx(i)) = corr(D(:,idx(i)),P(:,idx(i)));
        end
        s(isnan(s)) = 0;
        
        if strcmp(stat,'R2')
          s = s.^2;
        end
        
      case 'contingency'
        
        % return contingency table
        s = contingency(D,P);
        
      case 'confusion'
        
        % return the table of classification rates, i.e.,
        % the so-called confusion matrix
        
        tmp = contingency(D,P);
        tmp = tmp ./ repmat(sum(tmp,2),[1 size(tmp,2)]);
        tmp(isnan(tmp)) = 0;
        s = tmp;
        
      case 'binomial'
        
        k = sum(diag(contingency(D,P)));
        s = 1 - binocdf(k-1,nsamples,mean(D==mode(D)));
        
      case 'MAD'
        
        psi = abs(D-P);
        s =360*atan2(mean(sin(psi)),mean(cos(psi)))/(2*pi);
        
      case 'RMS'
        
        s = sqrt(mean((D-P).^2));
        
      case 'tlin'
        
        n=length(P);
        s=4*(sum(cos(P).*cos(D))*sum(sin(P).*sin(D))-sum(cos(P).*sin(D))*sum(sin(P).*cos(D)))/sqrt((n.^2-sum(cos(2*P)).^2-sum(sin(2*P)).^2)*(n.^2-sum(cos(2*D)).^2-sum(sin(2*D)).^2));
        
      case 'identity'
        
        % compute number of correctly identified trials by inspecting the
        % correlation between actual and predicted output
        
        [a,b] = max(corr(P',D'));
        s = mean(b==1:size(P,1));
        
      case 'expvar'
        
        s = (var(D) - var(D-P)) ./ var(D);
        
      otherwise
        error('unknown statistic %s',stat);
        
    end
    
    if minus, s = -s; end

  end
  
  
  function y = contingency(D,P)
    % compute contingency table with true class as rows and predicted class as columns
    
    % number of classes
    numcls = size(P,2);
    
    % predicted classes with the maximum posterior probability
    [tmp,pcls] = max(P,[],2);
    
    y = zeros(numcls);
    for i=1:numcls % true class is in rows
      ii = (D==i);
      for j=1:numcls % assigned class is in columns
        jj = (pcls==j);
        y(i,j) = sum(ii(:) & jj(:));
      end
    end
    
  end
end