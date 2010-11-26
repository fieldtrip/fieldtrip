classdef ft_mv_pca < ft_mv_preprocessor
%PCANALYZER performs a principal component analysis. Input data should be
%zero mean.
%
%   Options:
%    'proportion' : proportion of pc's or number of pc's. If < 1 then
%                 interpreted as a proportion of accounted variance; 
%                 otherwise as an absolute number; if empty 
%                 then all components are used (default = 0.80);
%
%   PARAMETERS:
%    pc     % principal components as column vectors
%    ev     % eigenvalues of principal components
%    accvar % cumulative variance accounted for per component
%
%   SEE ALSO:
%    princomp.m
%    pca2.m
%
%   REQUIRES:
%    statistics toolbox
%    pca toolbox
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
   
    proportion = []; % proportion of variance accounted for (0.80 is a good starting point)
  
    pc
    ev
    accvar
    
  end
  
  methods
    
    function obj = ft_mv_pca(varargin)
      
      obj = obj@ft_mv_preprocessor(varargin{:});

    end
    
    function obj = train(obj,X,Y)
      
      if nargin<3, Y = []; end
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      % data should be zero mean
      assert(all(abs(mean(X))<1e-6));
        
      if obj.verbose
        fprintf('estimating principal components\n');
      end
      
      if obj.proportion >= 1
        % use specialized fast approximation
        
        if obj.verbose
          fprintf('fast selection of %d principal components\n',obj.proportion);
        end

        % pca nameclashes with murphy toolbox and is renamed to pca2
        [U,S,V] = pca2(X,obj.proportion);
        obj.pc = V;
        
        if obj.verbose
          fprintf('quality of the approximation (estimate of spectral norm): %f\n',diffsnorm(X,U,S,V));
        end
                
      else
        
        % check availability
        if ~license('test','statistics_toolbox')
          error('requires Matlab statistics toolbox');
        end

        % in terms of variance accounted for
        
        [obj.pc,score,obj.ev] = princomp(X);
        obj.ev = obj.ev';
        
        % proportion of the variance that is accounted for
        obj.accvar = obj.ev/sum(obj.ev);
        
        % determine how many principal components to use
        if ~isempty(obj.proportion)
          
          if obj.proportion >= 1
            prop = 1:obj.proportion;
          else
            prop = 1:find(cumsum(obj.accvar) > obj.proportion,1,'first');
          end
          
          if obj.verbose
            fprintf('selected %d out of %d principal components\n',length(prop),size(obj.pc,2));
          end
          
          obj.pc = obj.pc(:,prop);
          obj.ev = obj.ev(prop);
          obj.accvar = obj.accvar(prop);
          
        end
        
      end
      
    end
    
    function Y = test(obj,X)
      
      Y = X * obj.pc;
      
    end
    
    function X = invert(obj,Y)
      
       X = Y * obj.pc';
     
    end
     
  end
end