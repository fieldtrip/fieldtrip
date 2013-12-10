classdef permutation
% PERMUTATION permutation testing class.
%
%   DESCRIPTION
%   Performs a permutation test on a dataset given a particular
%   crossvalidator object, a test statistic and a multivariate analysis.
%
%   EXAMPLE
%   m = dml.permutation('stat','accuracy','validator',dml.crossvalidator('mva',{dml.standardizer dml.naive}),'nperm',20,'verbose',true);
%   m = m.train(X,Y);
%   p = m.statistic
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  
  properties
    
    validator % crossvalidator object
    
    nperm = 100 % number of permutations
    
    outcome % outcome of the test statistic
    
    verbose = false; % whether or not to generate diagnostic output
    
    stat = ''; % statistic used to compute performance
        
  end
  
  methods
    
    function obj = permutation(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        else
          error('unrecognized fieldname %s',varargin{i});
        end
      end
      
      if isempty(obj.validator), error('specify cross-validator'); end
      
    end
    
    function obj = train(obj,X,Y)

      if iscell(Y)
        ndata = length(Y);
      else
        ndata = 1;
      end

      obj.outcome = zeros(1,obj.nperm+1);
      for i=1:obj.nperm
        
        % ensure random permutation with backward compatibility
        try
          RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
        catch
          rand('twister',sum(100*clock));
        end
        
        % create permuted data
        if ndata == 1
          Z = Y(randperm(size(Y,1)),:);
        else
          Z = cell(size(Y));
          for c=1:ndata
            Z{c} = Y{c}(randperm(size(Y{c},1)),:);
          end
        end
        
        if obj.verbose
          fprintf('testing permutation %d of %d: ',i,obj.nperm);
        end
        
        tmp = obj.validator.train(X,Z);
        obj.outcome(i) = tmp.statistic(obj.stat);
        if obj.verbose
          fprintf('%.2f\n',obj.outcome(i));
        end
        
      end
      obj.validator = obj.validator.train(X,Y);
      obj.outcome(end) = obj.validator.statistic(obj.stat);
      if obj.verbose
        fprintf('actual outcome: %.2f\n',obj.outcome(end));
      end
    end
    
    function p = statistic(obj)
      % return permutation test statistic; 
      
      T = obj.outcome(end); % real outcome
      S = sort(obj.outcome(1:(end-1))); % sampled outcomes
      
      N = length(obj.outcome);
      M = sum(S > T);
      
      % compute proportion of permutation distribution greater than or
      % equal to real statistic T
      m0 = sum(S == T);
      if m0
        p = (M + m0/2)/N; % deal with values equal to T
      else
        p = (M + 1)/N;
      end
            
    end
    
  end
  
end