classdef bootstrap
% BOOTSTRAP bootstrapping to determine parameter relevance.
%
%   DESCRIPTION
%   Performs a bootstrap test on a dataset given a particular
%   multivariate analysis. The stored results are the outputs of the
%   model function.
%
%   EXAMPLE
%   m = dml.bootstrap('mva',{dml.standardizer dml.naive},'nboot',100)
%   m = m.train(X,Y);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)
  
  properties
    
    mva % multivariate analysis
    
    nboot = 100 % number of permutations
    
    outcome; % bootstrap results
    
    verbose = false; % whether or not to generate diagnostic output
            
  end
  
  methods
    
    function obj = bootstrap(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        else
          error('unrecognized fieldname %s',varargin{i});
        end
      end
      
      if isempty(obj.mva), error('specify multivariate analysis'); end
      
      if ~isa(obj.mva,'dml.analysis'), obj.mva = dml.analysis(obj.mva); end
 
    end
    
    function obj = train(obj,X,Y)

      if iscell(Y)
        ndata = length(Y);
      else
        ndata = 1;
      end

      obj.outcome = cell(1,obj.nboot);
      for i=1:obj.nboot
        
        % ensure random permutation
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
        
        % create data by sampling with replacement
        if ndata == 1
          idx = randi(size(Y,1),1,size(Y,1));
          U = X(idx,:);
          V = Y(idx,:);
        else
          Z = cell(size(Y));
          for c=1:ndata
            idx = randi(size(Y{c},1),1,size(Y{c},1));
            U{c} = X{c}(idx,:);
            V{c} = Y{c}(idx,:);
          end
        end
        
        if obj.verbose
          fprintf('training bootstrap sample %d of %d\n',i,obj.nboot);
        end
        
        tmp = obj.mva.train(U,V);
        obj.outcome{i} = tmp.model;
        
      end
      
    end
    
    function [mu,se] = statistic(obj)
      % return mean and standard error of the mean for all parameters in
      % the model struct
      
      mu = obj.outcome{1}; se = mu;
      FN = fieldnames(mu);
      for f=1:length(FN)
        x = zeros([obj.nboot size(mu.(FN{f}))]); 
        for j=1:length(obj.outcome)
          x(j,:) = obj.outcome{j}.(FN{f});
        end
        mu.(FN{f}) = squeeze(mean(x));
        se.(FN{f}) = squeeze(std(x)/sqrt(obj.nboot));
      end
      
    end
    
  end
  
end