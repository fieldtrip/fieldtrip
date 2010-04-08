classdef pcanalyzer < preprocessor
%PCANALYZER performs a principal component analysis
%
%   Options:
%    'proportion' : proportion of pc's or number of pc's. If < 1 then
%                 interpreted as a proportion of accounted variance; 
%                 otherwise as an absolute number; if empty 
%                 then all components are used (default = 0.80);
%
%   NOTE:
%    pcanalyzer assumes the data has been standardized!
%
%   PARAMETERS:
%    accvar; % cumulative variance accounted for per component
%    pc; % principal components as column vectors
%    ev; % eigenvalues of principal components
%
%   SEE ALSO:
%    princomp.m
%
%   REQUIRES:
%    statistics toolbox
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
   proportion = []; % proportion of variance accounted for (0.80 is a good starting point)
  end
  
  methods
    
    function obj = pcanalyzer(varargin)
      
      % check availability
      if ~license('test','statistics_toolbox')
        error('requires Matlab statistics toolbox');
      end
      
      obj = obj@preprocessor(varargin{:});
    end
    
    function Y = map(obj,X)
      
      %M = dataset((obj.params.pc' * U.X')');
      Y = X * obj.params.pc;
      
    end
    
    function X = unmap(obj,Y)
      % CHECK!
      
      X = Y * obj.params.pc';
      
    end
    
    function p = estimate(obj,X,Y)
      
      if obj.proportion >= 1
        % use specialized fast approximation
        
        if obj.verbose
          fprintf('fast selection of %d principal components\n',obj.proportion);
        end

        [U,S,V] = pca(X,obj.proportion);
        p.pc = V;
                
      else
        % in terms of variance accounted for
        
        [p.pc,score,p.ev] = princomp(X);
        p.ev = p.ev';
        
        % proportion of the variance that is accounted for
        p.accvar = p.ev/sum(p.ev);
        
        % determine how many principal components to use
        if ~isempty(obj.proportion)
          
          if obj.proportion >= 1
            prop = 1:obj.proportion;
          else
            prop = 1:find(cumsum(p.accvar) > obj.proportion,1,'first');
          end
          
          if obj.verbose
            fprintf('selected %d out of %d principal components\n',length(prop),size(p.pc,2));
          end
          
          p.pc = p.pc(:,prop);
          p.ev = p.ev(prop);
          p.accvar = p.accvar(prop);
          
        end
        
      end
      
    end
    
  end
end